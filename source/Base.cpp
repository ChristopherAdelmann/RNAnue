#include "Base.hpp"

Base::Base(po::variables_map _params, std::string _subcall) : params(_params),
                                                              data(_params)
{

    std::string subcall = _params["subcall"].as<std::string>();

    // Define a vector of pairs from strings to member function pointers
    std::vector<std::pair<std::string, void (Data::*)()>> subcallMap = {
        {"preproc", &Data::preproc},
        {"align", &Data::align},
        {"detect", &Data::splitReadCalling},
        {"clustering", &Data::clustering},
        {"analysis", &Data::analysis}};

    // Check if the subcall is valid
    auto it = std::find_if(subcallMap.begin(), subcallMap.end(), [&](const auto &pair)
                           { return pair.first == subcall; });
    if (it == subcallMap.end() && subcall != "complete")
    {
        std::cout << "subcall: " << subcall << " invalid!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Call the appropriate function(s)
    for (auto &pair : subcallMap)
    {
        if (subcall == pair.first || subcall == "complete")
        {
            std::cout << "### Calling " << pair.first << "###" << std::endl;
            if (subcall == "complete")
            {
                std::cout << "### With complete workflow ###" << std::endl;
                data.setSubcall(pair.first);
            }
            (data.*(pair.second))();
        }
    }
}