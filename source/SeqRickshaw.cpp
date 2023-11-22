#include "SeqRickshaw.hpp"

SeqRickshaw::SeqRickshaw(const po::variables_map &params)
    : params(params)
{
    std::cout << helper::getTime() << " Started preprocessing" << std::endl;

    // retrieve parameters for the preprocessing
    minLen = params["minlen"].as<int>();
    phred = params["quality"].as<int>();
    wsize = params["wsize"].as<int>();

    cumDuration = std::chrono::duration<double>::zero();
    sumReads = 0;
}

void SeqRickshaw::setupLookupTables()
{
    std::cout << helper::getTime() << " Creating lookup table(s) for adapter trimming." << std::endl;

    // which adapters to trim (5'/3'/both)
    modus = params["modetrm"].as<int>();

    // readtype (SE or PE)
    readtype = params["readtype"].as<std::string>();

    const std::string adpt5File = params["adpt5"].as<std::string>();
    const std::string adpt3File = params["adpt3"].as<std::string>();

    // create folder for outputs of Rickshaw (e.g., lookup tables)
    fs::path outDirSub{params["outdir"].as<std::string>()};
    outDirSub /= fs::path(params["subcall"].as<std::string>());

    try
    {
        switch (modus)
        {
        case 0:
            adpt5Table = calcLookupTable("5\'-adapter", adpt5File);
            break;
        case 1:
        {
            adpt3Table = calcLookupTable("3\'-adapter", adpt3File);
            break;
        }
        case 2:
            adpt5Table = calcLookupTable("5\'-adapter", adpt5File);
            adpt3Table = calcLookupTable("3\'-adapter", adpt3File);
        }
    }
    catch (const seqan3::file_open_error &err)
    {
        std::cout << err.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (const std::exception &e)
    {
        std::cout << e.what() << std::endl;
    }

    if (params["savelookup"].as<bool>())
    {
        if (adpt3Table.has_value())
        {
            std::cout << helper::getTime() << " Saving lookup table for 3\'-adapter trimming." << std::endl;
            fs::path outDirSub3Adpt = outDirSub / fs::path("lookup_table_3adpts.txt");
            std::ofstream lookup3AdptOut(outDirSub3Adpt.string(), std::ios::trunc);
            writeLookupTables(lookup3AdptOut, adpt3Table.value());
        }

        if (adpt5Table.has_value())
        {
            std::cout << helper::getTime() << " Saving lookup table for 5\'-adapter trimming." << std::endl;
            fs::path outDirSub5Adpt = outDirSub / fs::path("lookup_table_5adpts.txt");
            std::ofstream lookup5AdptOut(outDirSub5Adpt.string(), std::ios::trunc);
            writeLookupTables(lookup5AdptOut, adpt5Table.value());
        }
    }

    std::cout << helper::getTime() << " Finished creating lookup table(s)" << std::endl;
}

// calculate the smart state transition table for the patterns
std::map<std::pair<std::string, std::string>, LookupTable> SeqRickshaw::calcLookupTable(std::string type, std::string adapterFile)
{
    std::cout << helper::getTime() << " Using: " << adapterFile << std::endl;

    Adapters tables;
    seqan3::sequence_file_input adaptersInput{adapterFile}; // read the adapters file

    for (auto &rec : adaptersInput)
    {
        // seqan3::debug_stream << "ID:  " << rec.id() << '\n';
        std::string readID = rec.id();
        auto readSeq = rec.sequence() | seqan3::views::to_char;

        LookupTable lookupTable = calcShift(readSeq);
        tables.emplace(std::make_pair(readID, std::string(readSeq.begin(), readSeq.end())), lookupTable);
    }
    return tables;
}

LookupTable SeqRickshaw::calcShift(auto &sequence)
{
    typedef std::tuple<long, long, long, std::bitset<1>> Entry;

    LookupTable lookup;

    // determine alphabet of sequence
    std::vector<char> alphabet = determineAlphabet(sequence);

    // intitial state - reading position to the right (e.g., 000000*)
    // create states set and add intial state
    std::string state(sequence.size(), '0');
    States states;

    std::size_t left = std::string::npos; // number of chars in the left block
    std::pair<std::size_t, std::size_t> right = std::make_pair(std::string::npos, std::string::npos);
    size_t readPos = sequence.size() - 1;

    states.push_back(std::make_tuple(state, left, right, readPos)); // left, right maintain info of previously matched chars

    int shift = 0;
    int match = 0; // match has been found

    std::string nextState;
    int nextStateID = 0;
    int nextReadPosVar = -1;

    int matched = 0; // how many
    size_t nrMatches = 0;

    int occ;

    std::string suffix;
    int suffixPos;

    char mismatch = ' ';

    std::string pattern;
    std::string subpat;

    size_t matchedCharsCount = 0;

    int goodSuffix = 0;
    int badChar = 0;

    States::size_type statesSize = states.size();
    for (States::size_type i = 0; i < statesSize; ++i)
    {
        state = std::get<0>(states[i]); // state (e.g., 000X000)
        left = std::get<1>(states[i]);  //
        right = std::get<2>(states[i]);
        readPos = std::get<3>(states[i]);

        for (unsigned j = 0; j < alphabet.size(); ++j)
        { // iterate through any other possibility

            std::size_t nextLeft = left;
            std::pair<std::size_t, std::size_t> nextRight = right;
            nextState = state;
            std::size_t nextReadPos;

            if (sequence[readPos] == alphabet[j])
            { //  match

                if (right.first == std::string::npos && right.second == std::string::npos)
                {
                    nextRight.first = sequence.size() - 1;
                    nextRight.second = sequence.size() - 1;
                }
                else
                {
                    if (readPos > right.second)
                    {
                        ++nextRight.second;
                    }
                    else
                    {
                        --nextRight.first;
                    }
                }
                nextState[readPos] = 'X';

                nextReadPos = calcReadPos(sequence, nextLeft, nextRight); // determines the readPos

                // exclude state with all matches
                if (std::count(nextState.begin(), nextState.end(), 'X') != sequence.size())
                {

                    nextStateID = addState(states, std::make_tuple(nextState, nextLeft, nextRight, nextReadPos), statesSize);

                    // add to lookup table
                    lookup.insert(std::make_pair(std::make_pair(i, alphabet[j]), std::make_tuple(0, nextStateID, nextReadPos, 0)));
                }
                else
                { // a complete match could have been detected
                    // add to lookup table
                    lookup.insert(std::make_pair(std::pair(i, alphabet[j]), std::make_tuple(0, 0, 0, 1)));
                }
            }
            else
            {
                // convert range to container
                auto tmp = sequence | ranges::to<std::vector<char>>;
                std::string pattern(tmp.begin(), tmp.end());

                // determine the suffix to check against the pattern
                if (right.first == std::string::npos && right.second == std::string::npos)
                {
                    suffix = alphabet[j]; // suffix consists only of mismatch
                }
                else
                {
                    suffix = pattern.substr(right.first, (right.second - right.first) + 1);
                    if (readPos < right.first)
                    {
                        suffix = alphabet[j] + suffix;
                    }
                    else
                    {
                        suffix = suffix + alphabet[j];
                    }
                }

                shift = transition(pattern, suffix, readPos, nextLeft, nextRight);

                nextState = std::string(sequence.size(), '0');

                // change stateID
                if (nextRight.first != std::string::npos && nextRight.second != std::string::npos)
                {
                    for (unsigned s = nextRight.first; s <= nextRight.second; ++s)
                    {
                        nextState[s] = 'X';
                    }
                }

                if (nextLeft != std::string::npos)
                {
                    for (unsigned l = 0; l <= nextLeft; ++l)
                    {
                        nextState[l] = 'X';
                    }
                }

                nextReadPos = calcReadPos(sequence, nextLeft, nextRight); // determines the readPos

                nextStateID = addState(states, std::make_tuple(nextState, nextLeft, nextRight, nextReadPos), statesSize);

                lookup.insert(std::make_pair(std::make_pair(i, alphabet[j]), std::make_tuple(shift, nextStateID, nextReadPos, 0)));
            }
        }
    }

    return lookup;
}

std::size_t SeqRickshaw::calcReadPos(auto &sequence, std::size_t &left, std::pair<std::size_t, std::size_t> &right)
{
    std::size_t readPos;

    //    std::cout << "readPos: ";

    // determine readPos
    if (right.first == std::string::npos && right.second == std::string::npos)
    {
        readPos = sequence.size() - 1; // start with first element
                                       //       std::cout << readPos << " - start with first element (to the right)" << std::endl;
    }
    else
    {
        // check if all positions are matches
        if ((right.second) - right.first == sequence.size() - 1 || left + (right.second - right.first) == sequence.size() - 2)
        {
            //          std::cout << "all matched " << std::endl;
            return std::string::npos;
            //            continue;
        }
        else
        {
            if (right.second + 1 < sequence.size())
            {                               // right end has not been reached
                readPos = right.second + 1; // go to the right
                                            //             std::cout << readPos << " - go to the right (there is still space) " << std::endl;
            }
            else
            { // continue to the left
                if (left == std::string::npos || left < right.first)
                {
                    readPos = right.first - 1;
                    //                std::cout << readPos << " - go to the left " << std::endl;
                }
            }
        }
    }

    return readPos;
}

int SeqRickshaw::transition(std::string pattern, std::string suffix, int readPos, std::size_t &left, std::pair<std::size_t, std::size_t> &right)
{
    int shift = 0;
    std::string subsuffix;

    std::size_t localShift = 0;

    // determine left most position of suffix + mismatch
    int leftMostPos = right.first;
    if (right.first == std::string::npos && right.second == std::string::npos)
    {
        leftMostPos = pattern.size() - 1;
    }
    else
    {
        if (readPos < right.first)
        {
            leftMostPos = right.first - 1;
        }
    }

    std::vector<std::size_t> occs;
    findAllOcc(occs, pattern, suffix);

    std::size_t fnd = std::string::npos;

    //
    for (unsigned i = 0; i < occs.size(); ++i)
    {
        localShift = leftMostPos - occs[i];
        if (left == std::string::npos || localShift - 1 >= left)
        {
            fnd = occs[i];
            break;
        }
    }

    // determine shift and the blocks (left, right)
    if (fnd != std::string::npos)
    {
        if (fnd != 0)
        {
            left = std::string::npos;
            right.first = fnd;
            right.second = fnd + (suffix.size() - 1);
            shift = leftMostPos - fnd;
        }
        else
        {
            left = suffix.size() - 1;
            right.first = std::string::npos;
            right.second = std::string::npos;
            shift = leftMostPos - fnd;
        }
    }
    else
    {
        for (unsigned j = 1; j < suffix.size(); ++j)
        {
            subsuffix = suffix.substr(j, suffix.size() - j);
            //            std::cout << "\t\tsubsuffix: " << subsuffix << std::endl;

            if (pattern.substr(0, subsuffix.size()) == subsuffix)
            {
                fnd = 0;
                //                std::cout << "\t\tfound " << std::endl;

                left = subsuffix.size() - 1;
                right.first = std::string::npos;
                right.second = std::string::npos;
                shift = leftMostPos + j;
            }
        }
        if (fnd == std::string::npos)
        {
            left = std::string::npos;
            right.first = std::string::npos;
            right.second = std::string::npos;
            shift = suffix.size();
        }
    }

    return shift;
}

/*
 * return the position within vector
 * */
int SeqRickshaw::addState(States &states, State state, States::size_type &size)
{
    // only add state if not in states
    int pos = -1;
    States::iterator it = std::find(states.begin(), states.end(), state);
    if (it != states.end())
    { // state is already within states
        pos = std::distance(states.begin(), it);
    }
    else
    {                            // state is not (yet) in states
        states.push_back(state); // push state to vector
        // stateIDs starting from 0
        pos = ++size - 1; //
    }
    return pos;
}

// determine next read position in state
int SeqRickshaw::nextReadPos(std::string state, int currReadPos)
{
    int found = -1;
    // check if readPos is not on the first
    if (currReadPos != state.size() - 1)
    {
        // search to the right of the current read position
        for (unsigned j = currReadPos + 1; j < state.size(); ++j)
        {
            if (state[j] == '0')
            { // found first
                std::cout << "\tfound in right at position: " << found << std::endl;
                found = j;
                return found;
            }
        }
    }

    // continue searching in the left part
    if (currReadPos != 0)
    {
        for (unsigned k = currReadPos - 1; k >= 0; --k)
        {
            if (state[k] == '0')
            {
                found = k;
                std::cout << "\tfound in left at position: " << found << std::endl;
                return found;
            }
        }
    }

    return -1;
}

/*
 * TODO: use action unique to implement this function more conventiently in c++20
 * */
std::vector<char> SeqRickshaw::determineAlphabet(auto _seq)
{
    std::vector<char> alphabet{};
    for (unsigned i = 0; i < _seq.size(); ++i)
    {
        if (std::find(alphabet.begin(), alphabet.end(), _seq[i]) == alphabet.end())
        {
            alphabet.push_back(_seq[i]);
        }
    }
    return alphabet;
}

void SeqRickshaw::findAllOcc(std::vector<std::size_t> &fnd, std::string str, std::string substr)
{
    std::size_t pos = str.find(substr);
    while (pos != std::string::npos)
    {
        fnd.push_back(pos);

        pos = str.find(substr, pos + substr.size());
    }
    std::reverse(fnd.begin(), fnd.end());
}

// the boyer moore search algorithm
std::size_t SeqRickshaw::boyermoore(auto &read, LookupTable tab, int m)
{
    std::size_t align = 0;
    int state = 0;
    int readPos = m - 1; // length of pattern to match against
    std::size_t shift = 0;

    std::tuple<int, int, int, int> entry;

    while (align < read.size() - m)
    {
        if (align + readPos >= read.size())
        {
            break;
        }
        auto c = (read | seqan3::views::to_char)[align + readPos];

        if (tab.find(std::make_pair(state, c)) == tab.end())
        {
            shift = 1;
            state = 0;
            readPos = m - 1;
        }
        else
        {
            entry = tab[std::make_pair(state, c)];
            shift = std::get<0>(entry);
            state = std::get<1>(entry);
            readPos = std::get<2>(entry);

            if (std::get<3>(entry) == 1)
            {
                return align;
            }
        }
        align = align + shift;
    }
    return std::string::npos;
}

//
std::pair<std::size_t, std::size_t> SeqRickshaw::trimming(auto &fwd)
{

    //  contains the search positons
    std::size_t fndPos;
    std::vector<size_t> res3Adpt;
    std::vector<size_t> res5Adpt;

    // marks the boundaries
    std::pair<std::size_t, std::size_t> bnds = std::make_pair(0, fwd.size());

    //  trim 5' adapters
    if (adpt5Table.has_value())
    {
        for (auto const &x : adpt5Table.value())
        {
            fndPos = boyermoore(fwd, x.second, x.first.second.size());

            if (fndPos == 0)
            { // adapter needs to be detected beginning with the first
                bnds.first = x.first.second.size();
            }
        }
    }

    // trim 3' adapters
    if (adpt3Table.has_value())
    {
        for (auto const &x : adpt3Table.value())
        {
            fndPos = boyermoore(fwd, x.second, x.first.second.size());
            if (fndPos < bnds.second)
            {
                bnds.second = fndPos;
            }
        }
    }
    return bnds;
}

SeqRickshaw::SeqRickshaw()
{
}

void SeqRickshaw::start(pt::ptree sample)
{
    setupLookupTables();

    pt::ptree input = sample.get_child("input");
    pt::ptree output = sample.get_child("output");

    std::string fwdRecInPath = input.get<std::string>("forward");
    std::string sampleName = std::filesystem::path(fwdRecInPath).stem().string();
    seqan3::sequence_file_input forwardRecIn{fwdRecInPath};

    std::cout << helper::getTime() << " Started processing " << sampleName << std::endl;

    std::pair<std::size_t, std::size_t> bndsFwd; //
    std::pair<std::size_t, std::size_t> bndsRev; //

    //
    if (readtype == "SE")
    {
        std::string outreads = output.get<std::string>("forward");

        std::ofstream myfile;
        myfile.open(outreads);

        for (auto &[seq, id, qual] : forwardRecIn)
        {

            bndsFwd = trimming(seq);
            // perform window trimming if specified
            if (params["wtrim"].as<std::bitset<1>>() == std::bitset<1>("1"))
            {
                // perform window trimming
                bndsFwd.second = nibble(seq, qual, bndsFwd);
            }

            auto trmReadFwd = seq | seqan3::views::slice(bndsFwd.first, bndsFwd.second);
            auto trmReadFwdQual = qual | seqan3::views::slice(bndsFwd.first, bndsFwd.second);

            // filter reads based on size
            if (std::ranges::size(trmReadFwd) != 0 && std::ranges::size(trmReadFwd) >= minLen)
            {
                auto bla = trmReadFwdQual | std::views::transform([](auto quality)
                                                                  { return seqan3::to_phred(quality); });
                auto sum = std::accumulate(bla.begin(), bla.end(), 0);
                auto qualscore = sum / std::ranges::size(bla);

                // std::cout << qualscore << std::endl;
                if (qualscore >= phred)
                {
                    myfile << "@" << id << '\n';
                    for (auto &s : trmReadFwd)
                    {
                        myfile << s.to_char();
                    }
                    myfile << '\n';
                    myfile << '+';
                    myfile << '\n';
                    for (auto &t : trmReadFwdQual)
                    {
                        myfile << t.to_char();
                    }
                    myfile << '\n';
                }
            }
        }
        myfile.close();
    }
    else
    { // readtype == "PE"
        std::string revRecInPath = input.get<std::string>("reverse");
        seqan3::sequence_file_input reverseRecIn{revRecInPath};

        // create output files for r1only/r2only - reads that don't survive the trimming

        std::string snglFwdOutPath = output.get<std::string>("R1only");
        std::string snglRevOutPath = output.get<std::string>("R2only");
        seqan3::sequence_file_output snglFwdOut{snglFwdOutPath};
        seqan3::sequence_file_output snglRevOut{snglRevOutPath};

        std::string mergeOutPath = output.get<std::string>("forward");
        std::ofstream mergeOut;
        mergeOut.open(mergeOutPath);

        int minOverlap = params["minovlps"].as<int>();

        auto start = std::chrono::high_resolution_clock::now();
        int reads = 0;
        Logger::log(LogLevel::INFO, sampleName + " ...");

        for (auto &&[rec1, rec2] : seqan3::views::zip(forwardRecIn, reverseRecIn))
        {
            if (reads % 1000 == 0)
            {
                auto end = std::chrono::high_resolution_clock::now();
                auto duration_per_read = std::chrono::duration_cast<std::chrono::duration<double>>(end - start) / reads;
                std::cout << helper::getTime() << " Reads processed: " << reads << " | Reads per second: " << 1 / duration_per_read.count() << std::endl;
            }

            trim3PolyG(rec1);
            trim3PolyG(rec2);
            const seqan3::dna5_vector adapterForward3 = "AGATCGGAAGAGC"_dna5;
            const seqan3::dna5_vector adapterReverse3 = "AGATCGGAAGAGC"_dna5;
            trimAdapter(adapterForward3, rec1, TrimConfig::Mode::THREE_PRIME);
            trimAdapter(adapterReverse3, rec2, TrimConfig::Mode::THREE_PRIME);

            std::string trmReadFwdID = rec1.id();
            std::string trmReadRevID = rec2.id();

            bool filtFwd = passesFilters(rec1);
            bool filtRev = passesFilters(rec2);

            if (filtFwd && filtRev)
            {
                auto mergedRecord = mergeRecordPair(rec1, rec2, minOverlap);
            }

            reads++;

            // Logger::log(LogLevel::INFO, rec1);
            // Logger::log(LogLevel::INFO, rec2);

            // if (filtFwd && filtRev)
            // {

            //     // std::pair<std::string, std::string> mrg = merging(trmFwdSeq, trmReadRev, trmFwdQual, trmReadRevQual);
            //     // if (mrg.first != "" && mrg.second != "")
            //     // {
            //     //     mergeOut << "@" << trmReadFwdID << '\n';
            //     //     mergeOut << mrg.first.c_str() << '\n';
            //     //     mergeOut << '+' << '\n';
            //     //     mergeOut << mrg.second << '\n';
            //     // }
            // }
            // else
            // {
            //     if (filtFwd)
            //     {
            //         // push to r1only
            //         // test_view | snglFwdOut;
            //         // snglFwdOut.emplace_back(trmFwdSeq, trmReadFwdID, trmFwdQual);
            //     }
            //     if (filtRev)
            //     {
            //         // push to r2only
            //         // snglRevOut.emplace_back(trmReadRev, trmReadRevID, trmReadRevQual);
            //     }
            // }
        }
        mergeOut.close();
    }
}

// window trimming method
std::size_t SeqRickshaw::nibble(auto &seq, auto &qual, std::pair<std::size_t, std::size_t> &bnds)
{
    std::size_t threePrimeEnd = bnds.second;

    while ((threePrimeEnd - bnds.first) >= wsize)
    {
        auto windowSeq = seq | seqan3::views::slice(threePrimeEnd - 3, threePrimeEnd);
        auto windowQual = qual | seqan3::views::slice(threePrimeEnd - 3, threePrimeEnd);

        // determine Phread score of window
        auto windowPhred = windowQual | std::views::transform([](auto quality)
                                                              { return seqan3::to_phred(quality); });
        auto windowPhredSum = std::accumulate(windowPhred.begin(), windowPhred.end(), 0);
        auto windowPhredScore = windowPhredSum / std::ranges::size(windowPhred);

        if (windowPhredScore >= phred)
        {
            break;
        }
        threePrimeEnd -= wsize;
    }
    return threePrimeEnd;
}

// determines the longest common substring between forward and reverse read
std::string SeqRickshaw::longestCommonSubstr(std::string s1, std::string s2)
{
    // Find length of both the strings.
    int m = s1.length();
    int n = s2.length();

    int dp[2][n + 1];
    int curr = 0, res = 0, end = 0;

    for (int i = 0; i <= m; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            if (i == 0 || j == 0)
            {
                dp[curr][j] = 0;
            }
            else
            {
                if (s1[i - 1] == s2[j - 1])
                {
                    dp[curr][j] = dp[1 - curr][j - 1] + 1;
                    if (res < dp[curr][j])
                    {
                        res = dp[curr][j];
                        end = i - 1;
                    }
                }
                else
                {
                    dp[curr][j] = 0;
                }
            }
        }
        curr = 1 - curr;
    }

    if (res == 0)
    {
        return "";
    }
    std::string ans;
    ans = s1.substr(end - res + 1, res);
    return ans;
}

/**
 * @brief Merges two records if they have a sufficient overlap.
 *
 * This function takes two records and checks if they have a sufficient overlap based on the specified minimum overlap value.
 * If the overlap is sufficient, the records are merged into a new record and returned as an optional. Otherwise, std::nullopt is returned.
 *
 * @tparam record_type The type of the records.
 * @param record1 The first record to be merged.
 * @param record2 The second record to be merged.
 * @param minOverlap The minimum required overlap between the records.
 * @return An optional containing the merged record if the overlap is sufficient, otherwise std::nullopt.
 */
template <typename record_type>
std::optional<record_type> SeqRickshaw::mergeRecordPair(record_type &record1, record_type &record2, int minOverlap)
{
    seqan3::align_cfg::method_global semiGlobalConfig{
        seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
        seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
        seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
        seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}};
    seqan3::align_cfg::scoring_scheme scoringSchemeConfig{
        seqan3::nucleotide_scoring_scheme{
            seqan3::match_score{1},
            seqan3::mismatch_score{-1}}};
    seqan3::align_cfg::gap_cost_affine gapSchemeConfig{
        seqan3::align_cfg::open_score{-2},
        seqan3::align_cfg::extension_score{-4}};

    auto outputConfig = seqan3::align_cfg::output_score{} |
                        seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_alignment{};

    auto alignmentConfig = semiGlobalConfig |
                           scoringSchemeConfig |
                           gapSchemeConfig |
                           outputConfig;

    auto &seq1 = record1.sequence();
    auto seq2ReverseComplement = record2.sequence() | std::views::reverse | seqan3::views::complement;

    for (auto const &result : seqan3::align_pairwise(std::tie(seq1, seq2ReverseComplement), alignmentConfig))
    {
        const int overlap = seq1.size() - result.sequence1_begin_position();
        if (overlap < minOverlap)
            continue;

        const int minScore = overlap - 2;
        if (result.score() >= minScore)
        {
            return constructMergedRecord(record1, record2, overlap);
        }
    }

    return std::nullopt;
}

/**
 * Constructs a merged record by combining two input records with a specified overlap.
 *
 * @tparam record_type The type of the input records.
 * @param record1 The first input record.
 * @param record2 The second input record.
 * @param overlap The length of the overlap between the two records.
 * @return The merged record.
 */
template <typename record_type>
record_type SeqRickshaw::constructMergedRecord(const record_type &record1, const record_type &record2, const int overlap)
{
    auto record2ReverseComplementSequence = record2.sequence() | seqan3::views::complement | std::views::reverse;
    auto record2ReverseComplementQualities = record2.base_qualities() | std::views::reverse;

    auto nonOverlapStartSequence = record1.sequence() | seqan3::views::slice(0, record1.sequence().size() - overlap);
    auto nonOverlapStartQualities = record1.base_qualities() | seqan3::views::slice(0, record1.sequence().size() - overlap);
    auto nonOverlapEndSequence = record2ReverseComplementSequence | seqan3::views::slice(overlap, record2ReverseComplementSequence.size());
    auto nonOverlapEndQualities = record2ReverseComplementQualities | seqan3::views::slice(overlap, record2ReverseComplementQualities.size());

    auto overlapSequenceRecord1 = record1.sequence() | seqan3::views::slice(record1.sequence().size() - overlap, record1.sequence().size());
    auto overlapSequenceRecord2 = record2ReverseComplementSequence | seqan3::views::slice(0, overlap);
    auto overlapQualitiesRecord1 = record1.base_qualities() | seqan3::views::slice(record1.sequence().size() - overlap, record1.sequence().size());
    auto overlapQualitiesRecord2 = record2ReverseComplementQualities | seqan3::views::slice(0, overlap);

    seqan3::dna5_vector overlapSequence;
    std::vector<seqan3::phred42> overlapQualities;

    overlapSequence.reserve(overlap);
    overlapQualities.reserve(overlap);

    for (auto const &[base1, base2, quality1, quality2] : seqan3::views::zip(overlapSequenceRecord1, overlapSequenceRecord2, overlapQualitiesRecord1, overlapQualitiesRecord2))
    {
        if (base1 == base2)
        {
            overlapSequence.push_back(base1);
            overlapQualities.push_back(std::max(quality1, quality2));
        }
        else
        {
            if (quality1 > quality2)
            {
                overlapSequence.push_back(base1);
                overlapQualities.push_back(quality1);
            }
            else
            {
                overlapSequence.push_back(base2);
                overlapQualities.push_back(quality2);
            }
        }
    }

    seqan3::dna5_vector mergedSequence;
    std::vector<seqan3::phred42> mergedQualities;

    mergedSequence.reserve(nonOverlapStartSequence.size() + overlapSequence.size() + nonOverlapEndSequence.size());
    mergedQualities.reserve(nonOverlapStartQualities.size() + overlapQualities.size() + nonOverlapEndQualities.size());

    mergedSequence.insert(mergedSequence.end(), std::make_move_iterator(nonOverlapStartSequence.begin()), std::make_move_iterator(nonOverlapStartSequence.end()));
    mergedSequence.insert(mergedSequence.end(), std::make_move_iterator(overlapSequence.begin()), std::make_move_iterator(overlapSequence.end()));
    mergedSequence.insert(mergedSequence.end(), std::make_move_iterator(nonOverlapEndSequence.begin()), std::make_move_iterator(nonOverlapEndSequence.end()));

    mergedQualities.insert(mergedQualities.end(), std::make_move_iterator(nonOverlapStartQualities.begin()), std::make_move_iterator(nonOverlapStartQualities.end()));
    mergedQualities.insert(mergedQualities.end(), std::make_move_iterator(overlapQualities.begin()), std::make_move_iterator(overlapQualities.end()));
    mergedQualities.insert(mergedQualities.end(), std::make_move_iterator(nonOverlapEndQualities.begin()), std::make_move_iterator(nonOverlapEndQualities.end()));

    record_type mergedRecord{std::move(mergedSequence), record1.id(), std::move(mergedQualities)};

    return mergedRecord;
}

std::pair<std::string, std::string> SeqRickshaw::merging(auto fwd, auto rev, auto fwdQual, auto revQual)
{
    auto forward = fwd | seqan3::views::to_char;
    auto reverse = rev | seqan3::views::complement | std::views::reverse | seqan3::views::to_char;

    auto forwardQual = fwdQual | seqan3::views::to_char;
    auto reverseQual = revQual | seqan3::views::to_char;

    std::string s1(forward.begin(), forward.end());
    std::string s2(reverse.begin(), reverse.end());

    std::string s1q(forwardQual.begin(), forwardQual.end());
    std::string s2q(reverseQual.begin(), reverseQual.end());

    std::string lcs = longestCommonSubstr(s1, s2);

    if (lcs.size() < 1 || lcs.size() < params["minovlps"].as<int>())
    {
        return std::make_pair("", "");
    }
    else
    {
        std::size_t s1Found = s1.find(lcs);
        std::size_t s2Found = s2.find(lcs);

        std::string mergedSeq = s1.substr(0, s1Found + lcs.size()) + s2.substr(0, s2Found);
        std::string mergedQual = s1q.substr(0, s1Found + lcs.size()) + s2q.substr(0, s2Found);

        return std::make_pair(mergedSeq, mergedQual);
    }
}

bool SeqRickshaw::passesFilters(auto &record)
{
    // Filter for mean quality
    auto phredQual = record.base_qualities() | std::views::transform([](auto q)
                                                                     { return q.to_phred(); });
    double sum = std::accumulate(phredQual.begin(), phredQual.end(), 0);
    double meanPhred = sum / std::ranges::size(phredQual);
    bool passesQual = meanPhred >= phred;

    // Filter for length
    bool passesLen = std::ranges::distance(record.sequence()) >= minLen;

    return passesQual && passesLen;
}

// write lookup table
void SeqRickshaw::writeLookupTables(std::ofstream &ofs, Adapters &adptLookTbls)
{
    for (auto const &[key, val] : adptLookTbls)
    {
        ofs << "ID: " << key.first << '\n';
        ofs << "Seq: " << key.second << '\n';
        for (auto const &[key2, val2] : val)
        {
            ofs << '(' << key2.first;
            ofs << ',' << key2.second << ") ->";
            ofs << '\t' << std::get<0>(val2);
            ofs << '\t' << std::get<1>(val2);
            ofs << '\t' << std::get<2>(val2);
            ofs << '\t' << std::get<3>(val2) << '\n';
        }
        ofs << '\n';
    }
}

void SeqRickshaw::processChunk(SeqRickshaw::FastqChunk const &chunk)
{
    // Your processing logic goes here.
    // Access the records using chunk.records.
    // Add any other processing steps.
    for (auto &&record : chunk.records)
    {
        // Wait to simulate processing.
        helper::simulateProcessing(std::chrono::microseconds(50));
    }
}

// Function to read a FASTQ file in chunks and process each chunk.
void SeqRickshaw::processFastqFileInChunks(std::string const &filename, size_t chunkSize, size_t numThreads)
{
    seqan3::sequence_file_input fastqFile{filename};

    // Mutex to synchronize access to shared data structures.
    std::mutex mutex;

    // Vector to store processed chunks for later use (optional).
    std::vector<FastqChunk> processedChunks;

    int chunkNumber = 0;

    // Function to read a chunk from the FASTQ file and process it.
    auto processChunkFunc = [&]()
    {
        while (true)
        {
            // Read a chunk of records from the FASTQ file.
            FastqChunk chunk;
            {
                std::lock_guard<std::mutex> lock(mutex);

                int count = 0;

                for (auto &&record : fastqFile)
                {
                    chunk.records.push_back(std::move(record));

                    if (++count == chunkSize)
                        break;
                }

                std::cout << "Chunk number: " << chunkNumber++ << '\n';
            }

            // Process the chunk.
            processChunk(chunk);

            // Optionally, store the processed chunk for later use.
            {
                std::lock_guard<std::mutex> lock(mutex);
                // Implement saving of preprocessed chunks here.
            }

            // Check if there are no more records in the file.
            if (fastqFile.begin() == fastqFile.end())
                break;
        }
    };

    // Create and run threads.
    std::vector<std::thread> threads;
    for (size_t i = 0; i < numThreads; ++i)
        threads.emplace_back(processChunkFunc);

    // Wait for all threads to finish.
    for (auto &thread : threads)
        thread.join();
}

/**
 * Trims adapter sequences from the given record.
 *
 * @param adapterSequence The adapter sequence to be trimmed.
 * @param record The record from which the adapter sequence will be trimmed.
 * @param trimmingMode The trimming mode to be applied.
 */
void SeqRickshaw::trimAdapter(const seqan3::dna5_vector &adapterSequence, auto &record, TrimConfig::Mode trimmingMode)
{
    seqan3::align_cfg::scoring_scheme scoringSchemeConfig{seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}};
    seqan3::align_cfg::gap_cost_affine gapSchemeConfig{seqan3::align_cfg::open_score{-1}, seqan3::align_cfg::extension_score{-2}};
    auto outputConfig = seqan3::align_cfg::output_score{} |
                        seqan3::align_cfg::output_begin_position{};

    auto alignment_config = TrimConfig::alignmentConfigFor(trimmingMode) |
                            scoringSchemeConfig |
                            gapSchemeConfig |
                            outputConfig;

    auto &seq = record.sequence();
    auto &qual = record.base_qualities();

    Logger::log(LogLevel::INFO, seq);

    for (auto const &result : seqan3::align_pairwise(std::tie(adapterSequence, seq), alignment_config))
    {
        int minScore = adapterSequence.size() - 2;
        if (result.score() >= minScore)
        {
            seq.erase(seq.begin() + result.sequence2_begin_position(), seq.end());
            qual.erase(qual.begin() + result.sequence2_begin_position(), qual.end());
        }
    }
}

/**
 * @brief Trims trailing polyG bases from the 3' end of a sequence record.
 *
 * This function removes any consecutive 'G' bases from the end of the sequence
 * if they meet the quality threshold. The mean quality threshold is set to a rank of 20
 * using the phred42 scale. If the number of consecutive polyG bases is greater than
 * or equal to 5, they are removed from both the sequence and the base qualities.
 *
 * @tparam record_type The type of the sequence record.
 * @param record The sequence record to trim.
 */
template <typename record_type>
void SeqRickshaw::trim3PolyG(record_type &record)
{
    auto &seq = record.sequence();
    auto &qual = record.base_qualities();

    std::vector<seqan3::phred42> qualitiesPhread;
    qualitiesPhread.reserve(qual.size());

    seqan3::phred42 qualityThreshold;
    qualityThreshold.assign_rank(30);

    auto seqIt = seq.rbegin();
    auto qualIt = qual.rbegin();

    auto sufficientMeanQuality = [&]()
    {
        auto qualities = qualitiesPhread | std::views::transform(
                                               [](auto quality)
                                               {
                                                   return seqan3::to_phred(quality);
                                               });

        auto sum = std::accumulate(qualities.begin(), qualities.end(), 0);
        return std::ranges::size(qualities) == 0 || sum / std::ranges::size(qualities) >= 20;
    };

    while (seqIt != seq.rend() && (*seqIt == 'G'_dna5 && sufficientMeanQuality()))
    {
        qualitiesPhread.push_back(*qualIt);
        ++seqIt;
        ++qualIt;
    }

    std::size_t polyGCount = qualitiesPhread.size();

    if (polyGCount >= 5)
    {
        seq.erase(seq.end() - polyGCount, seq.end());
        qual.erase(qual.end() - polyGCount, qual.end());
    }
}

/**
 * @brief Returns the semi-global alignment configuration for the given mode, assumes adapter as sequence1.
 *
 * @param mode The mode to use for trimming.
 * @return seqan3::align_cfg::method_global The semi-global alignment configuration.
 */
seqan3::align_cfg::method_global TrimConfig::alignmentConfigFor(TrimConfig::Mode mode)
{
    seqan3::align_cfg::method_global config;

    switch (mode)
    {
    case TrimConfig::Mode::FIVE_PRIME:
        config.free_end_gaps_sequence1_leading = true;
        config.free_end_gaps_sequence2_leading = true;
        config.free_end_gaps_sequence1_trailing = false;
        config.free_end_gaps_sequence2_trailing = true;
        break;
    case TrimConfig::Mode::THREE_PRIME:
        config.free_end_gaps_sequence1_leading = false;
        config.free_end_gaps_sequence2_leading = true;
        config.free_end_gaps_sequence1_trailing = true;
        config.free_end_gaps_sequence2_trailing = true;
        break;
    }

    return config;
}