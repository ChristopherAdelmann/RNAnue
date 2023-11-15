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
    typedef std::tuple<long,long,long,std::bitset<1>> Entry;

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

std::size_t SeqRickshaw::calcReadPos(auto& sequence, std::size_t& left, std::pair<std::size_t,std::size_t>& right) {
    std::size_t readPos;

//    std::cout << "readPos: "; 

    // determine readPos
    if(right.first == std::string::npos && right.second == std::string::npos) {
        readPos = sequence.size()-1; // start with first element
 //       std::cout << readPos << " - start with first element (to the right)" << std::endl;
    } else {
        // check if all positions are matches 
        if((right.second) - right.first == sequence.size() - 1 
                || left + (right.second - right.first) == sequence.size() - 2) {
  //          std::cout << "all matched " << std::endl;
            return std::string::npos;
//            continue;
        } else {
            if(right.second + 1 < sequence.size()) { // right end has not been reached
                readPos = right.second + 1; // go to the right
   //             std::cout << readPos << " - go to the right (there is still space) " << std::endl;
            } else { // continue to the left
                if(left == std::string::npos || left < right.first) {
                    readPos = right.first - 1;
    //                std::cout << readPos << " - go to the left " << std::endl;
                }
            } 
        }
    }

    return readPos;
}


        
int SeqRickshaw::transition(std::string pattern, std::string suffix, int readPos, std::size_t& left, std::pair<std::size_t,std::size_t>& right) {
    int shift = 0;
    std::string subsuffix;

    std::size_t localShift = 0;

    // determine left most position of suffix + mismatch
    int leftMostPos = right.first;
    if(right.first == std::string::npos && right.second == std::string::npos) {
        leftMostPos = pattern.size()-1;
    } else {
        if(readPos < right.first) {
            leftMostPos = right.first - 1;
        }
    }

    std::vector<std::size_t> occs;
    findAllOcc(occs, pattern, suffix);

    std::size_t fnd = std::string::npos;

    //
    for(unsigned i=0;i<occs.size();++i) {
        localShift = leftMostPos - occs[i];
        if(left == std::string::npos || localShift - 1 >= left) {
            fnd = occs[i];
            break;
        }
    }

    // determine shift and the blocks (left, right)
    if(fnd != std::string::npos) {
        if(fnd != 0) {
            left = std::string::npos;
            right.first = fnd;
            right.second = fnd + (suffix.size()-1);
            shift = leftMostPos - fnd;
        } else {
            left = suffix.size()-1;   
            right.first = std::string::npos;
            right.second = std::string::npos;
            shift = leftMostPos - fnd;
        }
    } else {
        for(unsigned j=1;j<suffix.size();++j) {
            subsuffix = suffix.substr(j,suffix.size()-j);
//            std::cout << "\t\tsubsuffix: " << subsuffix << std::endl;

            if(pattern.substr(0,subsuffix.size()) == subsuffix) {
                fnd = 0;
//                std::cout << "\t\tfound " << std::endl;

                left = subsuffix.size()-1;
                right.first = std::string::npos;
                right.second = std::string::npos;
                shift = leftMostPos + j;
            } 
        }
        if(fnd == std::string::npos) {
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
int SeqRickshaw::addState(States &states, State state, States::size_type &size) {
    // only add state if not in states
    int pos = -1;
    States::iterator it = std::find(states.begin(), states.end(), state);
    if(it != states.end()) { // state is already within states
        pos = std::distance(states.begin(), it);
    } else { // state is not (yet) in states
        states.push_back(state); // push state to vector
        // stateIDs starting from 0
        pos = ++size - 1; // 
    }
    return pos;
}


// determine next read position in state
int SeqRickshaw::nextReadPos(std::string state, int currReadPos) {
    int found = -1;
    // check if readPos is not on the first
    if(currReadPos != state.size()-1) { 
        // search to the right of the current read position
        for(unsigned j=currReadPos+1;j<state.size();++j) {
            if(state[j] == '0') { // found first 
                std::cout << "\tfound in right at position: " << found <<  std::endl;
                found = j;
                return found;
            }
        }
    }

    // continue searching in the left part
    if(currReadPos != 0) {
        for(unsigned k=currReadPos-1;k>=0;--k) {
            if(state[k] == '0') {
                found = k;
                std::cout << "\tfound in left at position: " <<  found << std::endl;
                return found;
            }
        }
    }

    return -1;
}

/*
 * TODO: use action unique to implement this function more conventiently in c++20
 * */
std::vector<char> SeqRickshaw::determineAlphabet(auto _seq) {
    std::vector<char> alphabet{};
    for(unsigned i=0; i < _seq.size(); ++i) {
        if(std::find(alphabet.begin(), alphabet.end(), _seq[i]) == alphabet.end()) {
            alphabet.push_back(_seq[i]);
        }
    }
    return alphabet;
}

void SeqRickshaw::findAllOcc(std::vector<std::size_t>& fnd, std::string str, std::string substr) {
    std::size_t pos = str.find(substr);
    while(pos != std::string::npos) {
        fnd.push_back(pos);

        pos = str.find(substr, pos + substr.size());
    }
    std::reverse(fnd.begin(),fnd.end());
}

// the boyer moore search algorithm
std::size_t SeqRickshaw::boyermoore(auto& read, LookupTable tab, int m) {
    std::size_t align = 0; 
    int state = 0;
    int readPos = m - 1; // length of pattern to match against
    std::size_t shift = 0;

    std::tuple<int,int,int,int> entry;

    while(align < read.size() - m) {
        if (align + readPos >= read.size())
        {
            break;
        }
        auto c = (read | seqan3::views::to_char)[align + readPos];

        if(tab.find(std::make_pair(state,c)) == tab.end()) {
            shift = 1;
            state = 0;
            readPos = m - 1;

        } else {
            entry = tab[std::make_pair(state,c)];
            shift = std::get<0>(entry);
            state = std::get<1>(entry); 
            readPos = std::get<2>(entry);
        
            if(std::get<3>(entry) == 1) {
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
    std::pair<std::size_t,std::size_t> bnds = std::make_pair(0,fwd.size());

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
            if(fndPos < bnds.second) {
                bnds.second = fndPos;
            }
        }
    }
    return bnds;
}

SeqRickshaw::SeqRickshaw() {
}

void SeqRickshaw::start(pt::ptree sample) {
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
    if(readtype == "SE") {
        std::string outreads = output.get<std::string>("forward");

        std::ofstream myfile;
        myfile.open(outreads);

        for (auto &[seq, id, qual] : forwardRecIn)
        {

            bndsFwd = trimming(seq);
            // perform window trimming if specified
            if(params["wtrim"].as<std::bitset<1>>() == std::bitset<1>("1")) {
			    // perform window trimming
			    bndsFwd.second = nibble(seq, qual, bndsFwd);
            }

            auto trmReadFwd = seq | seqan3::views::slice(bndsFwd.first,bndsFwd.second);
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
                    for(auto & s: trmReadFwd) {
						myfile << s.to_char();
					}
					myfile << '\n';
					myfile << '+';
					myfile << '\n';
					for(auto & t: trmReadFwdQual) {
						myfile << t.to_char();
					}
					myfile << '\n';
                }
            }
        }
        myfile.close();

    } else { // readtype == "PE"
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

        for (auto &&[rec1, rec2] : seqan3::views::zip(forwardRecIn, reverseRecIn))
        {
            bndsFwd = trimming(rec1.sequence());
            std::string trmReadFwdID = rec1.id();
            auto trmReadFwd = rec1.sequence() | seqan3::views::slice(bndsFwd.first, bndsFwd.second);
            auto trmReadFwdQual = rec1.base_qualities() | seqan3::views::slice(bndsFwd.first, bndsFwd.second);
            bool filtFwd = passesFilters(rec1);

            bndsRev = trimming(rec2.sequence());
            std::string trmReadRevID = rec2.id();
            auto trmReadRev = rec2.sequence() | seqan3::views::slice(bndsRev.first, bndsRev.second);
            auto trmReadRevQual = rec2.base_qualities() | seqan3::views::slice(bndsRev.first, bndsRev.second);
            bool filtRev = passesFilters(rec2);

            if(filtFwd && filtRev) {
                std::pair<std::string, std::string> mrg = merging(trmReadFwd,trmReadRev,trmReadFwdQual,trmReadRevQual);
                if(mrg.first != "" && mrg.second != "") {
                    mergeOut << "@" << trmReadFwdID << '\n';
                    mergeOut << mrg.first.c_str() << '\n';
                    mergeOut << '+' << '\n';
                    mergeOut << mrg.second << '\n';
                }

            } else {
                if(filtFwd) {
                    // push to r1only
                    snglFwdOut.emplace_back(trmReadFwd, trmReadFwdID, trmReadFwdQual);
                }
                if(filtRev) {
                    // push to r2only
                    snglRevOut.emplace_back(trmReadRev, trmReadRevID, trmReadRevQual);
                }
            }
        }
        mergeOut.close();
    }
}


// window trimming method
std::size_t SeqRickshaw::nibble(auto &seq, auto &qual, std::pair<std::size_t,std::size_t> &bnds) {
	std::size_t threePrimeEnd = bnds.second;

	while((threePrimeEnd - bnds.first) >= wsize ) {
		auto windowSeq = seq | seqan3::views::slice(threePrimeEnd-3,threePrimeEnd);
		auto windowQual = qual | seqan3::views::slice(threePrimeEnd-3,threePrimeEnd);

		// determine Phread score of window
		auto windowPhred = windowQual | std::views::transform([] (auto quality) { return seqan3::to_phred(quality); });
		auto windowPhredSum = std::accumulate(windowPhred.begin(), windowPhred.end(), 0);
		auto windowPhredScore = windowPhredSum / std::ranges::size(windowPhred);
		
		if(windowPhredScore >= phred) {
			break;
		}
		threePrimeEnd -= wsize;
	}
	return threePrimeEnd;
}

// determines the longest common substring between forward and reverse read
std::string SeqRickshaw::longestCommonSubstr(std::string s1, std::string s2) {
    // Find length of both the strings. 
    int m = s1.length(); 
    int n = s2.length(); 

    int dp[2][n+1];
    int curr=0,res=0,end=0;
    
    for(int i=0;i<=m;++i) {
        for(int j=0;j<=n;++j) {
            if(i==0 || j ==0) {
                dp[curr][j]=0;
            } else {
                if(s1[i-1] == s2[j-1]) {
                    dp[curr][j]=dp[1-curr][j-1]+1;
                    if(res<dp[curr][j]) {
                        res=dp[curr][j];
                        end=i-1;
                    }
                } else {
                    dp[curr][j]=0;
                }
            }
        }
        curr=1-curr;
    }

    if(res==0) {
        return "";
    }
    std::string ans;
    ans = s1.substr(end-res+1,res);
    return ans;

}


std::pair<std::string,std::string> SeqRickshaw::merging(auto fwd, auto rev, auto fwdQual, auto revQual) {
    auto forward = fwd | seqan3::views::to_char;
    auto reverse = rev | seqan3::views::complement | std::views::reverse | seqan3::views::to_char;

    auto forwardQual = fwdQual | seqan3::views::to_char;
    auto reverseQual = revQual | seqan3::views::to_char;

    std::string s1(forward.begin(),forward.end());
    std::string s2(reverse.begin(),reverse.end());

    std::string s1q(forwardQual.begin(),forwardQual.end());
    std::string s2q(reverseQual.begin(),reverseQual.end());

    std::string lcs = longestCommonSubstr(s1,s2);

    if(lcs.size() < 1 || lcs.size() < params["minovlps"].as<int>()) {
        return std::make_pair("","");
    } else {
        std::size_t s1Found = s1.find(lcs);
        std::size_t s2Found = s2.find(lcs);

        std::string mergedSeq = s1.substr(0,s1Found+lcs.size())+s2.substr(0,s2Found);
        std::string mergedQual = s1q.substr(0,s1Found+lcs.size())+s2q.substr(0,s2Found);

        return std::make_pair(mergedSeq,mergedQual);
    }
}

bool SeqRickshaw::passesFilters(auto &record)
{
    // Filter for mean quality
    auto qual = record.base_qualities() | std::views::transform([](auto q)
                                                                { return q.to_phred(); });
    double sum = std::accumulate(qual.begin(), qual.end(), 0);
    double meanPhred = sum / std::ranges::size(qual);
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
        for(auto const& [key2, val2] : val) {
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
