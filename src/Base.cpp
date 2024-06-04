#include "Base.hpp"

Base::Base(po::variables_map _params) : params(_params), data(_params) {
    std::string subcall = _params["subcall"].as<std::string>();
    // Define a vector of pairs from strings to member function pointers
    std::vector<std::pair<std::string, void (Data::*)()>> subcallMap = {
        {pi::PREPROCESS, &Data::preproc},
        {pi::ALIGN, &Data::align},
        {pi::DETECT, &Data::splitReadCalling},
        {pi::CLUSTER, &Data::clustering},
        {pi::ANALYZE, &Data::analysis}};

    // Check if the subcall is valid
    auto it =
        std::ranges::find_if(subcallMap, [&](const auto &pair) { return pair.first == subcall; });
    if (it == subcallMap.end() && subcall != pi::COMPLETE) {
        std::cout << "subcall: " << subcall << " invalid!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Call the appropriate function(s)
    for (auto &pair : subcallMap) {
        if (subcall == pair.first || subcall == pi::COMPLETE) {
            // Check if the preproc step is skipped
            if (!params[pi::PREPROCESS].as<bool>() && pair.first == pi::PREPROCESS) {
                continue;
            }

            if (subcall == pi::COMPLETE) {
                data.setSubcall(pair.first);
            }
            (data.*(pair.second))();
        }
    }
}
