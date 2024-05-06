#include "Analysis.hpp"

//
Analysis::Analysis(po::variables_map _params) : params(_params) {
    std::string line;
    std::ifstream anno;
    anno.open(params["features"].as<std::string>());

    // parse annotations
    if (!anno.is_open()) {
        perror("Error open");
        exit(EXIT_FAILURE);
    }
    while (getline(anno, line)) {
        if (line[0] == '#') {
            continue;
        }

        std::vector<std::string> tokens;
        std::istringstream iss(line);
        std::string token;
        while (std::getline(iss, token, '\t')) {
            tokens.push_back(token);
        }

        // boundaries
        std::pair<int, int> bnds = std::make_pair(std::stoi(tokens[3]), std::stoi(tokens[4]));
        std::pair<std::pair<int, int>, std::string> con =
            std::make_pair(bnds, "strand=" + tokens[6] + ";" + tokens[8]);

        auto it = features.find(tokens[0]);

        if (it == features.end()) {
            features.insert(
                std::pair<std::string, std::vector<std::pair<std::pair<int, int>, std::string>>>(
                    tokens[0], {con}));
        } else {
            features[tokens[0]].push_back(con);
        }
    }
}

//
void Analysis::createCountTable() {
    std::map<std::tuple<std::string, std::string, std::string, std::string>,
             std::vector<std::tuple<int, std::vector<float>, std::vector<float>>>>
        counts;

    std::ifstream intsfile;
    for (unsigned i = 0; i < interPaths.size(); ++i) {
        intsfile.open(interPaths[i]);

        if (!intsfile.is_open()) {
            perror("Error open");
            exit(EXIT_FAILURE);
        }

        std::string line;
        while (getline(intsfile, line)) {
            if (line[0] == '#') {
                continue;
            }

            std::vector<std::string> tokens;
            std::istringstream iss(line);
            std::string token;
            while (std::getline(iss, token, '\t')) {
                tokens.push_back(token);
            }

            // Create keys from tokens
            auto key = std::make_tuple(tokens[5], tokens[8], tokens[13], tokens[16]);

            // Check if key exists in counts
            auto it = counts.find(key);
            if (it == counts.end()) {
                // Key does not exist, create new entry
                std::vector<std::tuple<int, std::vector<float>, std::vector<float>>> content(
                    interPaths.size(), {0, {}, {}});

                std::get<0>(content[i]) = 1;
                std::get<1>(content[i]).push_back(std::stof(tokens[17]));
                std::get<2>(content[i]).push_back(std::stof(tokens[18]));

                counts.emplace(key, content);
            } else {
                // Key exists, update existing entry
                auto &oldVal = it->second[i];
                std::get<0>(oldVal) += 1;
                std::get<1>(oldVal).push_back(std::stof(tokens[17]));
                std::get<2>(oldVal).push_back(std::stof(tokens[18]));
            }
        }
    }

    // write back to file
    std::string outfile = params["outdir"].as<std::string>();
    fs::path outPath = fs::path(outfile) / "allints.txt";

    std::ofstream outFileHandle;
    outFileHandle.open(outPath.string());

    outFileHandle << "RNA1\tRNA2\tRNA1orientation\tRNA2orientation\t";

    for (const auto &interPath : interPaths) {
        std::string stem = fs::path(interPath).stem().string();
        outFileHandle << stem << "_counts\t";
        outFileHandle << stem << "_ges\t";
        outFileHandle << stem << "_gcs\t";
    }
    outFileHandle << std::endl;

    for (const auto &count : counts) {
        outFileHandle << std::get<0>(count.first) << "\t";
        outFileHandle << std::get<2>(count.first) << "\t";
        outFileHandle << std::get<1>(count.first) << "\t";
        outFileHandle << std::get<3>(count.first) << "\t";

        for (const auto &second : count.second) {
            outFileHandle << std::get<0>(second) << "\t";

            std::vector<float> vecNrg = std::get<1>(second);
            std::vector<float> vecCpl = std::get<2>(second);

            // determine ges
            float ges = 0.0;
            size_t vecNrgSize = vecNrg.size();
            if (vecNrgSize == 1) {
                ges = vecNrg[0];
            } else if (vecNrgSize > 1) {
                auto m = vecNrg.begin() + vecNrgSize / 2;
                std::nth_element(vecNrg.begin(), m, vecNrg.end());
                ges = vecNrg[vecNrgSize / 2];
            }
            outFileHandle << ges << "\t";

            // determine gcs
            float gcs = 0.0;
            size_t vecCplSize = vecCpl.size();
            if (vecCplSize == 1) {
                gcs = vecCpl[0];
            } else if (vecCplSize > 1) {
                auto m = vecCpl.begin() + vecCplSize / 2;
                std::nth_element(vecCpl.begin(), m, vecCpl.end());
                gcs = vecCpl[vecCplSize / 2];
            }
            outFileHandle << gcs << "\t";
        }
        outFileHandle << std::endl;
    }
    outFileHandle.close();
}

//
std::string Analysis::retrieveTagValue(std::string tags, std::string tagName,
                                       std::string oldValue) {
    std::size_t start_position = tags.find(tagName + "=");
    // gene name
    if (start_position != std::string::npos) {
        std::string sub = tags.substr(start_position + tagName.size() + 1, tags.length());
        std::size_t end_position = sub.find(";");
        oldValue = sub.substr(0, end_position);
    }
    return oldValue;
}

void Analysis::start(pt::ptree sample) {
    // retrieve input and output files
    pt::ptree input = sample.get_child("input");
    std::string splits = input.get<std::string>("splits");
    pt::ptree output = sample.get_child("output");
    std::string interactions = output.get<std::string>("interactions");

    interPaths.push_back(interactions);

    // input .sam record
    seqan3::sam_file_input fin{
        splits,
        seqan3::fields<seqan3::field::id, seqan3::field::flag, seqan3::field::ref_id,
                       seqan3::field::ref_offset, seqan3::field::seq, seqan3::field::tags>{}};

    std::vector<seqan3::sam_flag> flags;
    std::vector<std::string> refIDs;
    std::vector<std::optional<int32_t>> refOffsets;

    std::vector<size_t> ref_lengths{};
    for (auto &info : fin.header().ref_id_info) {
        ref_lengths.push_back(std::get<0>(info));
    }
    std::deque<std::string> ref_ids = fin.header().ref_ids();

    //    uint32_t flag, start, end;

    // open file
    std::ofstream outInts;
    outInts.open(interactions);

    std::string entry;  // stores output to write to file

    // variables for interactions file
    std::string qNAME, flag;
    uint32_t start, end;
    std::string geneName, product, annoStrand;
    float hybnrg, cmpl;

    seqan3::sam_tag_dictionary tags;

    outInts << "#QNAME\tSegment1Strand\tSegment1Start\tSegment1End\tSegment1RefName\tSegment1Name\t"
               "Segment1AnnoStrand\tSegment1Product\tSegment1Orientation\tSegment2Strand\tSegment2S"
               "tart\tSegment2End\tSegment2RefName\tSegment2Name\tSegment2AnnoStrand\tSegment2Produ"
               "ct\tSegment2Orientation\tenergy\tcomplementarity\n";

    int segCnt = 0;
    int segCntMatch = 0;
    int segFound = 0;

    int uniqueOverlaps = 0;
    int ambiguousOverlaps = 0;

    for (auto &&chimericReadPair : fin | seqan3::views::chunk(2)) {
        entry = "";
        hybnrg = 0.0;
        cmpl = 0.0;

        for (auto &read : chimericReadPair) {
            qNAME = read.id();
            flag =
                (static_cast<bool>(read.flag() & seqan3::sam_flag::on_reverse_strand)) ? "-" : "+";

            // start & end
            start = read.reference_position().value();
            end = start + read.sequence().size() - 1;

            // refID
            std::optional<int32_t> refIDidx = read.reference_id();
            std::string refID = ref_ids[refIDidx.value()];

            tags = read.tags();
            auto nrg = tags["XE"_tag];
            auto cpl = tags["XC"_tag];
            hybnrg = std::get<float>(nrg);
            cmpl = std::get<float>(cpl);

            if (features.count(refID) > 0) {
                std::vector<Feature> &featuresSubset = features[refID];

                geneName = ".";
                product = ".";
                annoStrand = ".";

                int overlaps = 0;
                // Rework this only keeps last feature not all features that could match
                for (unsigned i = 0; i < featuresSubset.size(); ++i) {
                    if ((start >= featuresSubset[i].first.first &&
                         start <= featuresSubset[i].first.second) ||
                        (end >= featuresSubset[i].first.first &&
                         end <= featuresSubset[i].first.second)) {
                        std::string oldGeneName = geneName;
                        geneName = retrieveTagValue(featuresSubset[i].second, "ID", geneName);
                        geneName = retrieveTagValue(featuresSubset[i].second, "gene", geneName);
                        product = retrieveTagValue(featuresSubset[i].second, "product", product);
                        annoStrand =
                            retrieveTagValue(featuresSubset[i].second, "strand", annoStrand);

                        if (oldGeneName != geneName) {
                            overlaps++;
                        }

                        segFound = 1;  // found a match (in annotation) for segment
                    }
                }

                if (overlaps > 1) {
                    Logger::log(LogLevel::WARNING, "Ambiguous annotation detected: ", refID, ":",
                                start, "-", end);
                    ambiguousOverlaps++;
                } else if (overlaps == 1) {
                    uniqueOverlaps++;
                }

                if (segFound == 1) {
                    if (segCnt == 0) {
                        entry += qNAME + "\t";
                    }
                    entry += flag + "\t";
                    entry += std::to_string(start) + "\t";
                    entry += std::to_string(end) + "\t";
                    entry += refID + "\t";
                    entry += geneName + "\t";
                    entry += annoStrand + "\t";
                    entry += "\"" + product + "\"\t";

                    if (flag == annoStrand) {
                        entry += "sense\t";
                    } else {
                        entry += "antisense\t";
                    }

                    segCntMatch++;
                    segFound = 0;
                }
            }
            segCnt++;  // on to the next segment
        }

        if (segCntMatch == 2) {
            outInts << entry << hybnrg << "\t" << cmpl << "\n";
        }

        segCnt = 0;  // reset segment
        segCntMatch = 0;
    }

    Logger::log(LogLevel::INFO, "Unique cluster annotations: ", uniqueOverlaps,
                "; Ambiguous cluster annotations: ", ambiguousOverlaps);

    outInts.close();
}
