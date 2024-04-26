#include "Align.hpp"

Align::Align(po::variables_map params) : params(params) {
  segemehlSysCall = params["segemehl"].as<std::string>();
  std::cout << "create align object" << std::endl;
  buildIndex();
}

void Align::alignReads(std::string query, std::string matched) {
  std::cout << "calling segemehl" << std::endl;

  std::string align = segemehlSysCall;
  align += " -S ";
  align += " -A " + std::to_string(params["accuracy"].as<int>());
  align += " -U " + std::to_string(params["minfragsco"].as<int>());
  align += " -W " + std::to_string(params["minsplicecov"].as<int>());
  align += " -Z " + std::to_string(params["minfraglen"].as<int>());
  align += " -t " + std::to_string(params["threads"].as<int>());
  align += " -m " + std::to_string(params["minlen"].as<int>());
  align += " -u /dev/null";
  align += " -r 1";
  align += " -i " + index;
  align += " -d " + params["dbref"].as<std::string>();
  align += " -q " + query;
  align += " -o " + matched;
  std::cout << align << std::endl;

  const char *call = align.c_str();
  system(call);
}

void Align::buildIndex() {
  // retrieve path of reference genome
  std::string ref = params["dbref"].as<std::string>();

  fs::path outDir = fs::path(params["outdir"].as<std::string>());
  fs::path gen = outDir / fs::path(ref).replace_extension(".idx").filename();

  if (fs::exists(gen)) {
    std::cout << "segemehl index found on filesystem\n";
  } else {
    std::cout << "generate index "
              << "\n";
    std::string genIndex = segemehlSysCall + " -x " + gen.string() + " -d " + ref;
    std::cout << genIndex << std::endl;
    const char *call = genIndex.c_str();
    system(call);
  }
  index = gen.string();
}

void Align::sortAlignments(std::string alignmentsPath) {
  Logger::log(LogLevel::INFO, "Sorting alignments");

  int const threads = params["threads"].as<int>();
  std::string tmpHeader = alignmentsPath + ".header";
  std::string tmpAlignments = alignmentsPath + ".alignments";

  // Command for extracting all header lines
  std::string headerCommand =
      "awk '/^@/ {print; next} {exit}' " + alignmentsPath + " > " + tmpHeader;
  // Command for sorting the alignments
  std::string sortCommand = "grep -v '^@' " + alignmentsPath + "| grep 'XJ:i:'" +
                            "| LC_ALL=C sort --parallel=" + std::to_string(threads) +
                            " -k1,1 -t$'\t' > " + tmpAlignments;
  // Command for merging the header and the sorted alignments
  std::string mergeCommand = "cat " + tmpHeader + " " + tmpAlignments + " > " + alignmentsPath;

  int result = system(headerCommand.c_str());
  if (result != 0) {
    Logger::log(LogLevel::ERROR, "Failed to execute header extraction");
    exit(1);
  }

  result = system(sortCommand.c_str());
  if (result != 0) {
    Logger::log(LogLevel::ERROR, "Failed to execute alignment sorting");
    exit(1);
  }

  result = system(mergeCommand.c_str());
  if (result != 0) {
    Logger::log(LogLevel::ERROR, "Failed to execute header and alignment merging");
    exit(1);
  }

  // Remove the temporary files
  fs::remove(tmpHeader);
  fs::remove(tmpAlignments);

  Logger::log(LogLevel::INFO, "Sorting alignments done");
}

//
seqan3::dna5 Align::string2dna5(std::string rna) {
  seqan3::dna5 seq{};
  for (auto &nt : rna) {
    seq.assign_char(nt);
  }
  return seq;
}

//
void Align::start(pt::ptree sample) {
  pt::ptree input = sample.get_child("input");
  pt::ptree output = sample.get_child("output");

  // align the reads
  std::string forward = input.get<std::string>("forward");
  std::string matched = output.get<std::string>("matched");
  alignReads(forward, matched);
  sortAlignments(matched);
}