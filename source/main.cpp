/***************************************************************************
*  Part of kASA: https://github.com/SilvioWeging/kASA
*
*  Copyright (C) 2020 Silvio Weging <silvio.weging@gmail.com>
*
*  Distributed under the Boost Software License, Version 1.0.
*  (See accompanying file LICENSE_1_0.txt or copy at
*  http://www.boost.org/LICENSE_1_0.txt)
**************************************************************************/
#include "kASA.hpp"
#include "modes/Read.hpp"
#include "modes/Shrink.hpp"
#include "modes/Compare.hpp"
#include "modes/Update.hpp"

#include "utils/dToStr.h"

//#include "Unittests.h"

namespace kASA_help {
	inline string getHelp(const string& m) {
		string out = "";
		if (m == "") {
			out = "Hello and welcome to kASA.\nYou did not specify any parameters.\nPlease read the README file on our Github https://github.com/SilvioWeging/kASA or select a mode and call <mode> -h or --help to get a list of parameters for the selected mode.\nPossible modes are: generateCF, build, identify, shrink, update, delete, getFrequency, redundancy, trie.\nGood luck and have a nice day!";
		}
		else if (m == "generateCF") {
			out = "This mode creates a content file out of genomic data with the help of the NCBI taxonomy.\n\
Necessary parameters:\n\
-i (--input) <file/folder>: Fasta file(s), can be gzipped. If you want to process multiple files at once, put them inside a folder and let the path end with `/`. No default.\n\
-u (--level) <level>: Taxonomic level at which you want to operate. All levels used in the NCBI taxonomy are available as well. To name a few: subspecies, species, genus, family, order, class, phylum, kingdom, superkingdom. Choose \"lowest\" if you want no linkage at a higher node in the taxonomic tree, this corresponds to other tools' \"sequence\" level. That means that no real taxid will be given and the name will be the line from the fasta containing the accession number.\n\
-f (--acc2tax) <folder or file>: As mentioned, either the folder containing the translation tables from accession number to taxid or a specific file. Can be gzipped.\n\
-y (--taxonomy) <folder>: This folder should contain the `nodes.dmp` and the `names.dmp` files.\n\
-c (--outgoing) <file>: Here, this parameter specifies where the content file should be written to.\n\
\n\
Optional parameters:\n\
-t (--temp) <path>: Path to temporary directory where files are stored that are either deleted automatically or can savely be deleted after kASA finishes. Defaults depend on your OS.\n\
-n (--threads) <number>: Number of parallel threads.Recommendation for different settings (due to IO bottleneck): HDD: 1, SSD: 2-4, RAM disk: 2-?. Note, that when compiling with C++17 enabled on Windows, some routines use all available cores provided by the hardware (because of the implementation of parallel STL algorithms). Default: 1.\n\
-m (--memory) <number>: Amount of Gigabytes available to kASA. If you don't provide enough, a warning will be written and it attempts to use as little as possible but may crash. If you provide more than your system can handle, it will crash or thrash. If you write \"inf\" instead of a number, kASA assumes that you have no memory limit. Default: 5 GB.\n\
-x (--callidx) <number>: Number given to this call of kASA so that no problems with temporary files occur if multiple instances of kASA are running at the same time. Default: 0.\n\
-v (--verbose): Prints out a little more information e.g.how much percent of your input was already read and analysed (if your input is not gzipped). Default: off.\n\
";
		}
		else if (m == "build") {
			out = "This mode creates a content file (if it doesn't already exist) and an index out of genomic data with the help of the NCBI taxonomy.\n\
Necessary parameters:\n\
-i (--input) <file/folder>: Fasta file(s), can be gzipped. If you want to process multiple files at once, put them inside a folder and let the path end with `/`. No default.\n\
-d (--database) <file>: Actually path and name of the index but let's call it database since `-i` was already taken...\n\
-c (--content) <file>: Path and name of the content file either downloaded or created from genomic data.\n\
If -c is not specified:\n\
-u (--level) <level>: Taxonomic level at which you want to operate. All levels used in the NCBI taxonomy are available as well. To name a few: subspecies, species, genus, family, order, class, phylum, kingdom, superkingdom. Choose \"lowest\" if you want no linkage at a higher node in the taxonomic tree, this corresponds to other tools' \"sequence\" level. That means that no real taxid will be given and the name will be the line from the fasta containing the accession number.\n\
-f (--acc2tax) <folder or file>: As mentioned, either the folder containing the translation tables from accession number to taxid or a specific file. Can be gzipped.\n\
-y (--taxonomy) <folder>: This folder should contain the `nodes.dmp` and the `names.dmp` files.\n\
\n\
Optional parameters:\n\
-a (--alphabet) <file> <number>: If you'd like to use a different translation alphabet formated in the NCBI compliant way, provide the file (gc.prt) and the id. Default: Hardcoded translation table.\n\
-z (--translated): Tell kASA, that the input consists of protein sequences. Currently in BETA.\n\
-t (--temp) <path>: Path to temporary directory where files are stored that are either deleted automatically or can savely be deleted after kASA finishes. Defaults depend on your OS.\n\
-n (--threads) <number>: Number of parallel threads.Recommendation for different settings (due to IO bottleneck): HDD: 1, SSD: 2-4, RAM disk: 2-?. Note, that when compiling with C++17 enabled on Windows, some routines use all available cores provided by the hardware (because of the implementation of parallel STL algorithms). Default: 1.\n\
-m (--memory) <number>: Amount of Gigabytes available to kASA.If you don't provide enough, a warning will be written and it attempts to use as little as possible but may crash. If you provide more than your system can handle, it will crash or thrash. If you write \"inf\" instead of a number, kASA assumes that you have no memory limit. Default: 5 GB.\n\
-x (--callidx) <number>: Number given to this call of kASA so that no problems with temporary files occur if multiple instances of kASA are running at the same time. Default: 0.\n\
-v (--verbose): Prints out a little more information e.g.how much percent of your input was already read and analysed (if your input is not gzipped). Default: off.\n\
--three: Use only three reading frames instead of six. Halves index size but may lead to artifacts due to additional reverse complement DNA inside some genomes. Default: off.\n\
";
		}
		else if (m == "identify") {
			out = "This mode analyzes genomic data for similarities with an index which was built with build.\n\
Necessary parameters:\n\
-i (--input) <file/folder>: Fastq or fasta file(s), can be gzipped. If you want to process multiple files at once, put them inside a folder and let the path end with `/`. No default.\n\
-d (--database) <file>: Actually path and name of the index but let's call it database since `-i` was already taken...\n\
-c (--content) <file>: Path and name of the content file either downloaded or created from genomic data.\n\
-p (--profile) <file>: Path and name of the profile that is put out.\n\
-q (--rtt) <file>: Path and name of the read ID to tax IDs output file. If not given, a profile-only version of kASA will be used which is much faster!\n\
\n\
Optional parameters:\n\
-r (--ram): Loads the index into primary memory. If you don't provide enough RAM for this, it will fall back to using secondary memory. Default: false.\n\
-k <upper> <lower>: Bounds for `k`, all `k`'s in between will be evaluated as well.If your intuition is more like `<lower > <upper>` then that's okay too. Default: 12 7.\n\
--kH <upper>: Set only the upper bound.\n\
--kL <lower>: Set only the lower bound.\n\
-b (--beasts) <number>: Number of hit taxa shown for each read. Default: 3.\n\
--json: Sets the output format to json. Default.\n\
--jsonl: Sets the output format to json lines.\n\
--tsv: Sets the output format to a tab separated, per line format.\n\
--kraken: Sets the output format to a kraken like tsv format.\n\
-e (--unique): Ignores duplicates of k-mers in every read.This helps removing bias from repeats but messes with error scores, so consider it BETA.\n\
-a (--alphabet) <file> <number>: If you'd like to use a different translation alphabet formated in the NCBI compliant way, provide the file (gc.prt) and the id. Default: Hardcoded translation table.\n\
-z (--translated): Tell kASA, that the input consists of protein sequences. Currently in BETA.\n\
\n\
-t (--temp) <path>: Path to temporary directory where files are stored that are either deleted automatically or can savely be deleted after kASA finishes. Defaults depend on your OS.\n\
-n (--threads) <number>: Number of parallel threads.Recommendation for different settings (due to IO bottleneck): HDD: 1, SSD: 2-4, RAM disk: 2-?. Note, that when compiling with C++17 enabled on Windows, some routines use all available cores provided by the hardware (because of the implementation of parallel STL algorithms). Default: 1.\n\
-m (--memory) <number>: Amount of Gigabytes available to kASA.If you don't provide enough, a warning will be written and it attempts to use as little as possible but may crash. If you provide more than your system can handle, it will crash or thrash. If you write \"inf\" instead of a number, kASA assumes that you have no memory limit. Default: 5 GB.\n\
-x (--callidx) <number>: Number given to this call of kASA so that no problems with temporary files occur if multiple instances of kASA are running at the same time. Default: 0.\n\
-v (--verbose): Prints out a little more information e.g.how much percent of your input was already read and analysed (if your input is not gzipped). Default: off.\n\
--threshold <float>: Set a minimum relative score so that everything below it will not be included in the output. Default: 0.0.\n\
--six: Use all six reading frames instead of three. Doubles number of input k-mers but avoids artifacts due to additional reverse complement DNA inside some genomes. Default: off.\n\
";
		}
		else if (m == "update") {
			out = "This mode updates an existing index with new genomic data.\n\
Necessary parameters:\n\
-i (--input) <file/folder>: Fasta file(s), can be gzipped. If you want to process multiple files at once, put them inside a folder and let the path end with `/`. No default.\n\
-d (--database) <file>: Current index. No default.\n\
-c (--content) <file>: Path and name of the content file either downloaded or created from genomic data.\n\
-u (--level) <level>: Taxonomic level at which you want to operate. All levels used in the NCBI taxonomy are available as well. To name a few: subspecies, species, genus, family, order, class, phylum, kingdom, superkingdom. Choose \"lowest\" if you want no linkage at a higher node in the taxonomic tree, this corresponds to other tools' \"sequence\" level. That means that no real taxid will be given and the name will be the line from the fasta containing the accession number.\n\
-f (--acc2tax) <folder or file>: As mentioned, either the folder containing the translation tables from accession number to taxid or a specific file. Can be gzipped.\n\
-y (--taxonomy) <folder>: This folder should contain the `nodes.dmp` and the `names.dmp` files.\n\
-o (--outgoing) <file>: Either the existing index or a new file name, depending on whether you want to keep the old file or not. Default: overwrite.\n\
\n\
Optional parameters:\n\
-a (--alphabet) <file> <number>: If you'd like to use a different translation alphabet formated in the NCBI compliant way, provide the file (gc.prt) and the id. Default: Hardcoded translation table.\n\
-z (--translated): Tell kASA, that the input consists of protein sequences. Currently in BETA.\n\
\n\
-t (--temp) <path>: Path to temporary directory where files are stored that are either deleted automatically or can savely be deleted after kASA finishes. Defaults depend on your OS.\n\
-n (--threads) <number>: Number of parallel threads.Recommendation for different settings (due to IO bottleneck): HDD: 1, SSD: 2-4, RAM disk: 2-?. Note, that when compiling with C++17 enabled on Windows, some routines use all available cores provided by the hardware (because of the implementation of parallel STL algorithms). Default: 1.\n\
-m (--memory) <number>: Amount of Gigabytes available to kASA.If you don't provide enough, a warning will be written and it attempts to use as little as possible but may crash. If you provide more than your system can handle, it will crash or thrash. If you write \"inf\" instead of a number, kASA assumes that you have no memory limit. Default: 5 GB.\n\
-x (--callidx) <number>: Number given to this call of kASA so that no problems with temporary files occur if multiple instances of kASA are running at the same time. Default: 0.\n\
-v (--verbose): Prints out a little more information e.g.how much percent of your input was already read and analysed (if your input is not gzipped). Default: off.\n\
--three: Use only three reading frames instead of six. Halves index size but may lead to artifacts due to additional reverse complement DNA inside some genomes. Default: off.\n\
";
		}
		else if (m == "shrink") {
		out = "This mode reduces the size of an existing index.\n\
Necessary parameters:\n\
-d (--database) <file>: Index file. No default.\n\
-c (--content) <file>: Path and name of the content file either downloaded or created from genomic data.\n\
-s (--strategy) <1, 2 or 3>: Shrink the index in the first, second or third way. Default is 2.\n\
-g (--percentage) <integer>: Deletes the given percentage of k-mers from every taxon.This parameter may also be applied when building the index (for example: -g 50 skips every second k-mer).\n\
-o (--outgoing) <file>: Output path and name of your shrunken index file. Your other index cannot be overwritten with this. Default: takes your index file and appends a \"_s\".\n\
\n\
Optional parameters:\n\
-t (--temp) <path>: Path to temporary directory where files are stored that are either deleted automatically or can savely be deleted after kASA finishes. Defaults depend on your OS.\n\
-n (--threads) <number>: Number of parallel threads.Recommendation for different settings (due to IO bottleneck): HDD: 1, SSD: 2-4, RAM disk: 2-?. Note, that when compiling with C++17 enabled on Windows, some routines use all available cores provided by the hardware (because of the implementation of parallel STL algorithms). Default: 1.\n\
-m (--memory) <number>: Amount of Gigabytes available to kASA.If you don't provide enough, a warning will be written and it attempts to use as little as possible but may crash. If you provide more than your system can handle, it will crash or thrash. If you write \"inf\" instead of a number, kASA assumes that you have no memory limit. Default: 5 GB.\n\
-x (--callidx) <number>: Number given to this call of kASA so that no problems with temporary files occur if multiple instances of kASA are running at the same time. Default: 0.\n\
-v (--verbose): Prints out a little more information e.g.how much percent of your input was already read and analysed (if your input is not gzipped). Default: off.\n\
";
		}
		else if (m == "delete") {
		out = "This mode deletes deprecated entries from an existing index via the delnodes.dmp file from the NCBI taxonomy.\n\
Necessary parameters:\n\
-d (--database) <file>: Index file. No default.\n\
-c (--content) <file>: Path and name of the content file either downloaded or created from genomic data.\n\
-l (--deleted) <file>: delete taxa via the NCBI taxonomy file.\n\
-o (--outgoing) <file>: Either the existing index or a new file name, depending on whether you want to keep the old file or not. Default: overwrite.\n\
Optional parameters:\n\
-t (--temp) <path>: Path to temporary directory where files are stored that are either deleted automatically or can savely be deleted after kASA finishes. Defaults depend on your OS.\n\
-n (--threads) <number>: Number of parallel threads.Recommendation for different settings (due to IO bottleneck): HDD: 1, SSD: 2-4, RAM disk: 2-?. Note, that when compiling with C++17 enabled on Windows, some routines use all available cores provided by the hardware (because of the implementation of parallel STL algorithms). Default: 1.\n\
-m (--memory) <number>: Amount of Gigabytes available to kASA.If you don't provide enough, a warning will be written and it attempts to use as little as possible but may crash. If you provide more than your system can handle, it will crash or thrash. If you write \"inf\" instead of a number, kASA assumes that you have no memory limit. Default: 5 GB.\n\
-x (--callidx) <number>: Number given to this call of kASA so that no problems with temporary files occur if multiple instances of kASA are running at the same time. Default: 0.\n\
-v (--verbose): Prints out a little more information e.g.how much percent of your input was already read and analysed (if your input is not gzipped). Default: off.\n\
";
		}
		else if (m == "getFrequency") {
		out = "This mode generates the frequencies file (containing the number of k-mers for every taxon) for an existing index.\n\
Necessary parameters:\n\
-d (--database) <file>: Index file. No default.\n\
-c (--content) <file>: Path and name of the content file either downloaded or created from genomic data.\n\
\n\
Optional parameters:\n\
-t (--temp) <path>: Path to temporary directory where files are stored that are either deleted automatically or can savely be deleted after kASA finishes. Defaults depend on your OS.\n\
-n (--threads) <number>: Number of parallel threads.Recommendation for different settings (due to IO bottleneck): HDD: 1, SSD: 2-4, RAM disk: 2-?. Note, that when compiling with C++17 enabled on Windows, some routines use all available cores provided by the hardware (because of the implementation of parallel STL algorithms). Default: 1.\n\
-m (--memory) <number>: Amount of Gigabytes available to kASA.If you don't provide enough, a warning will be written and it attempts to use as little as possible but may crash. If you provide more than your system can handle, it will crash or thrash. If you write \"inf\" instead of a number, kASA assumes that you have no memory limit. Default: 5 GB.\n\
-x (--callidx) <number>: Number given to this call of kASA so that no problems with temporary files occur if multiple instances of kASA are running at the same time. Default: 0.\n\
-v (--verbose): Prints out a little more information e.g.how much percent of your input was already read and analysed (if your input is not gzipped). Default: off.\n\
";
		}
		else if (m == "redundancy") {
		out = "This mode calculates for an existing index, how many taxa belong to each k-mer.\n\
Necessary parameters:\n\
-d (--database) <file>: Index file. No default.\n\
-c (--content) <file>: Path and name of the content file either downloaded or created from genomic data.\n\
\n\
Optional parameters:\n\
-t (--temp) <path>: Path to temporary directory where files are stored that are either deleted automatically or can savely be deleted after kASA finishes. Defaults depend on your OS.\n\
-x (--callidx) <number>: Number given to this call of kASA so that no problems with temporary files occur if multiple instances of kASA are running at the same time. Default: 0.\n\
";
		}
		else if (m == "trie") {
		out = "This mode creates the prefix trie file needed for \"identify\".\n\
Necessary parameters:\n\
-d (--database) <file>: Index file. No default.\n\
\n\
Optional parameters:\n\
-t (--temp) <path>: Path to temporary directory where files are stored that are either deleted automatically or can savely be deleted after kASA finishes. Defaults depend on your OS.\n\
-x (--callidx) <number>: Number given to this call of kASA so that no problems with temporary files occur if multiple instances of kASA are running at the same time. Default: 0.\n\
";
		}
		return out;
	}
}

int main(int argc, char* argv[]) {
	try {

		vector<string> vParameters(argv, argv + argc);

		string cMode = "", sDBPathOut = "", sTempPath = "", sInput = "", contentFileIn = "", readToTaxaFile = "", tableFile = "", indexFile = "", delnodesFile = "", codonTable = "", sTaxonomyPath = "", sAccToTaxFiles = "", sTaxLevel = "", sStxxlMode = "", sCodonID = "1";
		bool bSpaced = false, bVerbose = false, bTranslated = false, bRAM = false, bUnique = false, bUnfunny = false, bSixFrames = false, bThreeFrames = false;
		kASA::Shrink::ShrinkingStrategy eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::TrieHalf;
		kASA::Compare::OutputFormat eOutputFormat = kASA::Compare::OutputFormat::Json;
		int32_t iNumOfThreads = 1, iHigherK = 12, iLowerK = 7, iNumOfCall = 0, iNumOfBeasts = 3;
		uint64_t iMemorySizeAvail = 0;
		float fPercentageOfThrowAway = 0.f, threshold = 0.f;
		uint8_t iTrieDepth = 6, iPrefixCheckMode = 0;

		auto timeRightNow = chrono::system_clock::to_time_t(chrono::system_clock::now());
		cout << "OUT: " << "kASA version " << kASA_VERSION << " ran on " <<
#if  _WIN64
			"Windows "
#elif __linux__
			"Linux "
#elif __APPLE__
			"macOS "
#else
			"Something other than Win/Linux/macOS "
#endif
			<< "at " << ctime(&timeRightNow) << "OUT: ";
		for (int32_t i = 0; i < argc; ++i) {
			 cout << vParameters[i] << " ";
		}
		cout << endl;

		// Default temporary path:
#if  _WIN64
		sTempPath = string(getenv("appdata")) + "/";
#elif __linux__
		sTempPath = "/var/tmp/";
#elif __APPLE__
		sTempPath = "/private/tmp/";
#else
		cerr << "Default temporary path not found, please specify it if you didn't already..." << endl;
#endif

		if (argc == 1) {
			cout << "OUT: \n\n" << kASA_help::getHelp("") << endl;
			return 0;
			//throw runtime_error("No Parameters given!");
		}
		cMode = vParameters[1];

		// still available parameters: w

		for (int32_t i = 2; i < argc; ++i) {
			string sParameter = vParameters[i];
			if (sParameter == "-h" || sParameter == "--help") {
				cout << "OUT: HELP CALLED FOR MODE " << cMode << "\n\n" << kASA_help::getHelp(cMode) << endl;
				return 0;
			}
			else if (sParameter == "-o" || sParameter == "--outgoing") {
				sDBPathOut = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-t" || sParameter == "--temp") {
				sTempPath = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-u" || sParameter == "--level") {
				sTaxLevel = Utilities::removeSpaceAndEndline(vParameters[++i]);
				if (sTaxLevel == "sequence") {
					sTaxLevel = "lowest";
				}
			}
			else if (sParameter == "-e" || sParameter == "--unique") {
				bUnique = true;
			}
			else if (sParameter == "-f" || sParameter == "--acc2tax") {
				sAccToTaxFiles = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-y" || sParameter == "--taxonomy") {
				sTaxonomyPath = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-v" || sParameter == "--verbose") {
				bVerbose = true;
			}
			else if (sParameter == "-z" || sParameter == "--translated") {
				bTranslated = true;
			}
			else if (sParameter == "-j" || sParameter == "--sloppy") {
				bUnfunny = true;
#ifdef _DEBUG
				//kASA::kASA::setAAToAATable(Utilities::removeSpaceAndEndline(vParameters[++i]));
#endif
			}
			else if (sParameter == "-d" || sParameter == "--database") {
				indexFile = Utilities::removeSpaceAndEndline(vParameters[++i]);
				if (cMode != "build") {
					if (!ifstream(indexFile)) {
						throw runtime_error("Index file not found");
					}
					if (!ifstream(indexFile + "_info.txt")) {
						throw runtime_error("Info file not found");
					}
					if (!ifstream(indexFile + "_f.txt") && cMode == "identify") {
						throw runtime_error("Frequency file not found");
					}
					if (!ifstream(indexFile + "_trie") && cMode == "identify") {
						throw runtime_error("Trie file not found");
					}
					if (!ifstream(indexFile + "_trie.txt") && cMode == "identify") {
						throw runtime_error("Trie info file not found");
					}
				}
			}
			else if (sParameter == "-a" || sParameter == "--alphabet") {
				codonTable = Utilities::removeSpaceAndEndline(vParameters[++i]);
				sCodonID = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-b" || sParameter == "--beasts") {
				iNumOfBeasts = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				if (iNumOfBeasts == 0) {
					iNumOfBeasts = 1;
				}
			}
			else if (sParameter == "-r" || sParameter == "--ram") {
				bRAM = true;
			}
			else if (sParameter == "-g" || sParameter == "--percentage") {
				fPercentageOfThrowAway = stof(Utilities::removeSpaceAndEndline(vParameters[++i]));
			}
			else if (sParameter == "-x" || sParameter == "--callidx") {
				iNumOfCall = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
			}
			else if (sParameter == "-n" || sParameter == "--threads") {
				iNumOfThreads = abs(stoi(Utilities::removeSpaceAndEndline(vParameters[++i])));
			}
			else if (sParameter == "-k") {
				iHigherK = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				iLowerK = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				iHigherK = (iHigherK > 12) ? 12 : iHigherK;
				iLowerK = (iLowerK < 1) ? 1 : iLowerK;
				if (iLowerK > iHigherK) {
					swap(iLowerK, iHigherK);
				}
			}
			else if (sParameter == "--kH") {
				iHigherK = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				iHigherK = (iHigherK > 12) ? 12 : iHigherK;
			}
			else if (sParameter == "--kL") {
				iLowerK = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				iLowerK = (iLowerK < 1) ? 1 : iLowerK;
			}
			else if (sParameter == "-i" || sParameter == "--input") {
				sInput = Utilities::removeSpaceAndEndline(vParameters[++i]);
				if (!ifstream(sInput) && sInput.back() != '/') {
					throw runtime_error("Input file not found");
				}
			}
			else if (sParameter == "-q" || sParameter == "--rtt") {
				readToTaxaFile = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-p" || sParameter == "--profile") {
				tableFile = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-m" || sParameter == "--memory") {
				const auto& userGivenMemory = Utilities::removeSpaceAndEndline(vParameters[++i]);
				if (userGivenMemory == "inf") {
					iMemorySizeAvail = numeric_limits<uint64_t>::max() / (1024ull * 1024ull);
				}
				else {
					iMemorySizeAvail = 1024ull * stoi(userGivenMemory);
				}
			}
			else if (sParameter == "-s" || sParameter == "--strategy") {
				int32_t iChoice = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				switch (iChoice) {
				case 1:
					eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::EveryNth;
					break;
				case 2:
					eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::TrieHalf;
					break;
				case 3:
					eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::Entropy;
					break;
				case 4:
					eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::Overrepresented;
					break;
				default:
					eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::TrieHalf;
					break;
				}
			}
			else if (sParameter == "-c" || sParameter == "--content") {
				contentFileIn = Utilities::removeSpaceAndEndline(vParameters[++i]);
				if (!ifstream(contentFileIn) && cMode != "build" && cMode != "generateCF") {
					throw runtime_error("Content file not found");
				}
			}
			else if (sParameter == "-l" || sParameter == "--deleted") {
				delnodesFile = Utilities::removeSpaceAndEndline(vParameters[++i]);
				if (!ifstream(delnodesFile)) {
					throw runtime_error("Deleted nodes file not found");
				}
			}
			else if (sParameter == "--json") {
				eOutputFormat = kASA::Compare::OutputFormat::Json;
			}
			else if (sParameter == "--jsonl") {
				eOutputFormat = kASA::Compare::OutputFormat::JsonL;
			}
			else if (sParameter == "--tsv") {
				eOutputFormat = kASA::Compare::OutputFormat::tsv;
			}
			else if (sParameter == "--kraken") {
				eOutputFormat = kASA::Compare::OutputFormat::Kraken;
			}
			else if (sParameter == "--stxxl") {
				sStxxlMode = vParameters[++i];
			}
			else if (sParameter == "--array") {
				iPrefixCheckMode = 2;
			}
			else if (sParameter == "--trie") {
				iPrefixCheckMode = 1;
			}
			else if (sParameter == "--table") {
				iPrefixCheckMode = 3;
			}
			else if (sParameter == "--six") {
				bSixFrames = true;
			}
			else if (sParameter == "--three") {
				bThreeFrames = true;
			}
			else if (sParameter == "--threshold") {
				threshold = stof(vParameters[++i]);
			}
			else {
				throw runtime_error("Some unknown parameter has been inserted, please check your command line.");
			}
		}
#ifdef ENVIRONMENT32
		iMemorySizeAvail = (iMemorySizeAvail == 0 || iMemorySizeAvail >= numeric_limits<uint32_t>::max() || iMemorySizeAvail >= 2048) ? 1024 : iMemorySizeAvail; // TODO: wenn kleiner als 1024, auf 1024 setzen, sonst bug in stxxl::sort
#else
		iMemorySizeAvail = (iMemorySizeAvail == 0) ? 5120 : iMemorySizeAvail;
#endif
		iMemorySizeAvail *= 1024ull * 1024ull;

#if _WIN32 || _WIN64
		MEMORYSTATUSEX statex;
		statex.dwLength = sizeof(statex);
		GlobalMemoryStatusEx(&statex);
		if (statex.ullTotalPhys < iMemorySizeAvail) {
			cout << "OUT: WARNING! The requested memory of " << iMemorySizeAvail/(1024ull*1024ull*1024ull) << "GB RAM is too much for your system (" << statex.ullTotalPhys / (1024ull * 1024ull * 1024ull) <<"GB RAM in total). kASA may crash or slow down to a crawl..." << endl;
		}

#elif __GNUC__ || __clang__
		unsigned long long pages = sysconf(_SC_PHYS_PAGES);
		unsigned long long page_size = sysconf(_SC_PAGE_SIZE);
		unsigned long long iMaxPhysMemory = pages * page_size;

		if (iMaxPhysMemory < iMemorySizeAvail) {
			cout << "OUT: WARNING! The requested memory of " << iMemorySizeAvail / (1024ull * 1024ull * 1024ull) << "GB RAM is too much for your system (" << iMaxPhysMemory / (1024ull * 1024ull * 1024ull) << "GB RAM in total). kASA may crash or slow down to a crawl..." << endl;
		}
#endif

		if (thread::hardware_concurrency() < static_cast<uint32_t>(iNumOfThreads)) {
			cout << "OUT: WARNING! The requested number of " << iNumOfThreads << " threads is too much for your system (" << thread::hardware_concurrency() << " cores available). kASA may crash or slow down to a crawl..." << endl;
		}

		if (cMode == "build") {
			kASA::Read kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, bTranslated, sStxxlMode, !bThreeFrames, bUnfunny);
			if (codonTable != "") {
				kASAObj.setCodonTable(codonTable, sCodonID);
			}
			auto start = std::chrono::high_resolution_clock::now();

			// No content file yet created
			if (contentFileIn == "") {
				if (sTaxLevel != "lowest" && (sAccToTaxFiles == "" || sTaxonomyPath == "")) {
					throw runtime_error("No acc2Tax file or taxonomy path given...");
				}
				else {
					contentFileIn = indexFile + "_content.txt";
					if (bVerbose) {
						cout << "OUT: Creating content file: " << contentFileIn << endl;
					}
					kASAObj.generateContentFile(sTaxonomyPath, sAccToTaxFiles, sInput, contentFileIn, sTaxLevel);
				}
			}
			if (iMemorySizeAvail*0.9 < 1024ull * 1024ull * 1024ull) {
				throw runtime_error("Not enough memory given!");
			}
			kASAObj.BuildAll(contentFileIn, sInput, indexFile, static_cast<uint64_t>(iMemorySizeAvail*0.9 - 1024ull * 1024ull * 1024ull), fPercentageOfThrowAway);

			auto end = std::chrono::high_resolution_clock::now();
			cout << "OUT: Time: " << chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << endl;
		}
		else if (cMode == "generateCF") {
			if (contentFileIn == "") {
				throw runtime_error("Where should I put the content file?");
			}
			if (sTaxLevel != "lowest" && (sAccToTaxFiles == "" || sTaxonomyPath == "")) {
				throw runtime_error("No acc2Tax file or taxonomy path given...");
			}
			kASA::kASA kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose);
			kASAObj.generateContentFile(sTaxonomyPath, sAccToTaxFiles, sInput, contentFileIn, sTaxLevel);
		}
		else if (cMode == "update") {
			kASA::Update kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, bTranslated, sStxxlMode, !bThreeFrames);
			if (codonTable != "") {
				kASAObj.setCodonTable(codonTable, sCodonID);
			}

			// Content file was created together with the index
			if (contentFileIn == "") {
				contentFileIn = indexFile + "_content.txt";
			}

			// get taxIds that are new and map the accession numbers to those
			if (sAccToTaxFiles != "" && sTaxonomyPath != "") {
				kASAObj.addToContentFile(sTaxonomyPath, sAccToTaxFiles, sInput, sTaxLevel, contentFileIn);
			}

			auto start = std::chrono::high_resolution_clock::now();
			kASAObj.UpdateFromFasta(contentFileIn, indexFile, sInput, sDBPathOut, (indexFile == sDBPathOut) || (sDBPathOut == ""), static_cast<uint64_t>(iMemorySizeAvail*0.9 - 1024ull * 1024ull * 1024ull), fPercentageOfThrowAway);
			auto end = std::chrono::high_resolution_clock::now();
			cout << "OUT: Time: " << chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << endl;
		}
		else if (cMode == "delete") {

			if (sDBPathOut == "") {
				throw runtime_error("No output file given!");
			}

			kASA::Update kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, bTranslated, sStxxlMode);
			auto start = std::chrono::high_resolution_clock::now();
			kASAObj.DeleteFromLib(contentFileIn, indexFile, sDBPathOut, delnodesFile, (indexFile == sDBPathOut));
			auto end = std::chrono::high_resolution_clock::now();
			cout << "OUT: Time: " << chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << endl;
		}
		else if (cMode == "shrink") {
			if (indexFile == sDBPathOut) {
				throw runtime_error("Paths and names of input and output are the same!");
			}

			// Default value
			bool bCopyContentFile = false;
			if (sDBPathOut == "") {
				sDBPathOut = indexFile + "_s";
			}

			// Content file was created together with the index
			if (contentFileIn == "") {
				contentFileIn = indexFile + "_content.txt";
				bCopyContentFile = true;
			}

			kASA::Shrink kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, sStxxlMode);
			kASAObj.ShrinkLib(indexFile, sDBPathOut, eShrinkingStrategy, contentFileIn, fPercentageOfThrowAway);

			if (bCopyContentFile) {
				// Copy file so that both indices can be used by default parameters for -c
				ifstream CFile(contentFileIn, std::ios::binary);
				ofstream SCFile(sDBPathOut + "_content.txt", std::ios::binary);
				if (!SCFile || !CFile) {
					throw runtime_error("Content file couldn't be opened for reading/writing!");
				}
				SCFile << CFile.rdbuf();
			}
		}
		else if (cMode == "identify") {
			kASA::Compare kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, iNumOfBeasts, bVerbose, bTranslated, sStxxlMode, bSixFrames, bUnfunny);
			if (codonTable != "") {
				kASAObj.setCodonTable(codonTable, sCodonID);
			}

			// Content file was created together with the index
			if (contentFileIn == "") {
				contentFileIn = indexFile + "_content.txt";
			}

			kASAObj.format = eOutputFormat;
			auto start = std::chrono::high_resolution_clock::now();
			kASAObj.CompareWithLib_partialSort(contentFileIn, indexFile, sInput, readToTaxaFile, tableFile, iTrieDepth, iMemorySizeAvail, bSpaced, bRAM, bUnique, iPrefixCheckMode, threshold);
			auto end = std::chrono::high_resolution_clock::now();
			cout << "OUT: Time: " << chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << endl;
#if _WIN32 || _WIN64
			//system("PAUSE"); //DEBUG
#endif
		}
		else if (cMode == "getFrequency") {
			kASA::kASA kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, sStxxlMode);
			if (contentFileIn == "") {
				contentFileIn = indexFile + "_content.txt";
			}

			kASAObj.GetFrequencyK(contentFileIn, indexFile, indexFile + "_f.txt");
		}
		else if (cMode == "redundancy") {
			// Content file was created together with the index
			if (contentFileIn == "") {
				contentFileIn = indexFile + "_content.txt";
			}

			kASA::Shrink kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, sStxxlMode);
			uint32_t iIdxCounter = 1;
			ifstream content(contentFileIn);
			string sDummy = "";
			while (getline(content, sDummy)) {
				if (sDummy != "") {
					++iIdxCounter;
				}
			}

			ifstream fLibInfo(indexFile + "_info.txt");
			uint64_t iSizeOfLib = 0, shrunken = 0;
			fLibInfo >> iSizeOfLib;
			fLibInfo >> shrunken;
			fLibInfo.close();
			
			if (shrunken) {
				throw runtime_error("redundancy cannot be called on shrunken indices!");
			}

			stxxlFile libFile(indexFile, stxxlFile::RDONLY);
			unique_ptr<const contentVecType_32p> libVec(new const contentVecType_32p(&libFile, iSizeOfLib));
			const auto& iCutoffNumber = kASAObj.histogram(libVec,iIdxCounter);
			if (iCutoffNumber == 1) {
				cout << "OUT: 99% of the k-mers in your index have only one taxon. Using unique frequencies makes sense." << endl;
			}
			else {
				if (iCutoffNumber < 4) {
					cout << "OUT: 99% of the k-mers in your index have " << iCutoffNumber << " or less taxa. Using unique frequencies could make sense." << endl;
				}
				else {
					cout << "OUT: 99% of the k-mers in your index have " << iCutoffNumber << " or less taxa. You should consider looking at the non-unique frequencies as well." << endl;
				}
			}
		}
		else if (cMode == "trie") {
			kASA::kASA kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, sStxxlMode);

			ifstream fLibInfo(indexFile + "_info.txt");
			uint64_t iSizeOfLib = 0;
			fLibInfo >> iSizeOfLib;
			fLibInfo.close();

			auto start = std::chrono::high_resolution_clock::now();

			stxxlFile libFile(indexFile, stxxlFile::RDONLY);
			const contentVecType_32p libVec(&libFile, iSizeOfLib);

			//cout << stxxl::is_sorted(libVec.cbegin(), libVec.cend(), [](const packedBigPair& a, const packedBigPair& b) {return (a.first < b.first || (!(b.first < a.first) && a.second < b.second)); }) << endl;

			Trie T(static_cast<int8_t>(12), static_cast<int8_t>(iLowerK), iTrieDepth);
			T.SaveToStxxlVec(&libVec, indexFile);

			auto end = std::chrono::high_resolution_clock::now();
			cout << "OUT: Time: " << chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << endl;
		}
		else if (cMode == "half") {
			if (indexFile == sDBPathOut) {
				throw runtime_error("Paths and names of input and output are the same!");
			}
			kASA::Shrink kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, sStxxlMode);
			//kASAObj.GetFrequencyK(contentFileIn, databaseFile, sDBPathOut + "_f.txt");
			kASAObj.ShrinkLib(indexFile, sDBPathOut, kASA::Shrink::ShrinkingStrategy::TrieHalf, contentFileIn);
		}
		else if (cMode == "debug") {
			//kASA::Compare UnitTests(sTempPath, 1, 12, 9, 0, iNumOfBeasts);
			//UnitTests.testC2V_1();
			//UnitTests.testC2V_2();
			//UnitTests.testC2V_3();
			//UnitTests.testC2V_4();
			//UnitTests.testC2V_5();
			//UnitTests.testC2V_6();
			//UnitTests.testC2V_7();
			//UnitTests.testC2V_8();
			//UnitTests.testC2V_9();
		}
		else if (cMode == "test") {
			stxxlFile temp(indexFile, stxxl::file::RDONLY);
			ifstream fLibInfo(indexFile + "_info.txt");
			uint64_t iSizeOfLib = 0;
			fLibInfo >> iSizeOfLib;

			ifstream searchFile(vParameters[2]);
			string kMerString = "";
			vector<uint64_t> vSearchVec;
			while (getline(searchFile, kMerString)) {
				vSearchVec.push_back(kASA::kASA::aminoacidTokMer(kMerString));
			}
			uint32_t iSearchCounter = 0;

			const contentVecType_32p tempVec(&temp, iSizeOfLib);
			for (const auto& entry : tempVec) {
				const auto& kMer = entry.first;
				if (kMer < vSearchVec[iSearchCounter]) {
					continue;
				}
				else {
					if (!(vSearchVec[iSearchCounter] < kMer)) {
						cout << kASA::kASA::kMerToAminoacid(kMer, 12) << " " << entry.second << endl;
					}
					else {
						iSearchCounter++;
						if (iSearchCounter >= vSearchVec.size()) {
							break;
						}
						else {
							if (!(vSearchVec[iSearchCounter] < kMer)) {
								cout << kASA::kASA::kMerToAminoacid(kMer, 12) << " " << entry.second << endl;
							}
						}
					}
				}
			}
		}
		else if (cMode == "howmuchtaxids") {
			stxxlFile temp(indexFile, stxxl::file::RDONLY);
			ifstream fLibInfo(indexFile + "_info.txt");
			uint64_t iSizeOfLib = 0;
			fLibInfo >> iSizeOfLib;

			ofstream searchFile(sTempPath+"frequentkMers.txt");
			uint64_t iSeenkMer = 0;
			uint32_t iSearchCounter = 0;
			uint32_t iMagicNumber = 4;
			set<uint32_t> taxIDs;

			const contentVecType_32p tempVec(&temp, iSizeOfLib);
			for (const auto& entry : tempVec) {
				const auto& kMer = entry.first;
				if (kMer == iSeenkMer) {
					taxIDs.insert(entry.second);
					iSearchCounter++;
				}
				else {
					if (iSearchCounter >= iMagicNumber) {
						searchFile << kASA::kASA::kMerToAminoacid(kMer, 12);
						for (const auto& elem : taxIDs) {
							searchFile << " " << elem;
						}
						searchFile << endl;
					}
					iSearchCounter = 0;
					taxIDs.clear();
					iSeenkMer = kMer;
				}
			}
		}
		else if (cMode == "showVec") {
			stxxlFile temp(indexFile, stxxl::file::RDONLY);
			ifstream sizeFile(indexFile + "_info.txt");
			uint64_t iSize = 0;
			sizeFile >> iSize;
			bool index_t = false;
			sizeFile >> index_t;
			if (index_t) {
				kASA::kASA::showVec(index_t_p(&temp, iSize));
			}
			else {
				kASA::kASA::showVec(contentVecType_32p(&temp, iSize));
			}
		}
		else if (cMode == "transform") {
			stxxlFile temp(indexFile, stxxl::file::RDONLY);

			ifstream sizeFile(indexFile + "_info.txt"); // "_info.txt"
			uint64_t iSize = 0;
			sizeFile >> iSize;

			const contentVecType_32p tempV(&temp, iSize);

			Utilities::createFile(sDBPathOut);
			Utilities::createFile(sDBPathOut + "_2");

			stxxlFile tempOut(sDBPathOut, stxxl::file::RDWR);
			stxxl::VECTOR_GENERATOR<uint64_t, 4U, 4U, 2101248, stxxl::RC>::result tempPV(&tempOut, iSize);

			stxxlFile tempOut2(sDBPathOut+"_2", stxxl::file::RDWR);
			stxxl::VECTOR_GENERATOR<uint32_t, 4U, 4U, 2101248, stxxl::RC>::result tempPV2(&tempOut2, iSize);
			
			ofstream countsFile(sDBPathOut+"_counts.txt");

			auto t1It = tempV.cbegin();
			auto tOIt = tempPV.begin();
			auto tO2It = tempPV2.begin();
			uint64_t iCount = 0;
			uint64_t iSeen = 0;
			for (; t1It != tempV.cend();) {
				if (t1It->first == iSeen) {
					*tO2It = t1It->second;
					++t1It;
				}
				else {
					*tOIt = iSeen = t1It->first;
					countsFile << iCount << endl;
					*tO2It = t1It->second;
					++t1It; ++tOIt;
				}
				++tO2It;
				++iCount;
			}

			tempPV.resize(tOIt - tempPV.begin(), true);

			ofstream outSizeFile(sDBPathOut + "_info.txt");
			outSizeFile << tempPV.size() << endl << iSize;
			tempPV.export_files("_");
			tempPV2.export_files("_");
		}
		else if (cMode == "fuckit") {


			uint32_t iIdxCounter = 1;
			unordered_map<uint32_t, uint32_t> mIDsAsIdx; mIDsAsIdx[0] = 0;
			ifstream content(contentFileIn);
			string sDummy = "";
			while (getline(content, sDummy)) {
				if (sDummy != "") {
					const auto& line = Utilities::split(sDummy, '\t');
					if (line.size() == 4) {
						mIDsAsIdx[stoul(line[1])] = iIdxCounter;
						++iIdxCounter;
					}
				}
			}

			ifstream sizeFile(indexFile + "_info.txt"); // "_info.txt"
			uint64_t iSize = 0;
			sizeFile >> iSize;

			stxxlFile temp(indexFile, stxxl::file::RDONLY);
			const contentVecType_32p tempV(&temp, iSize);

			Utilities::createFile(sTempPath + "_tmp");

			stxxlFile tempOut(sTempPath + "_tmp", stxxl::file::RDWR);
			contentVecType_32p tempPV(&tempOut, iSize);


			auto t1It = tempV.cbegin();
			auto tOIt = tempPV.begin();
			while (t1It != tempV.cend()) {
				tOIt->second = t1It->second;
				uint64_t tVal = 0;
				for (int32_t i = 55, j = 0; i >= 5; i -= 10, j += 5) {
					tVal |= (t1It->first & (31ULL << i)) << j;
				}
				tOIt->first = tVal;

				++t1It;
				++tOIt;
			}

			struct SCompareStructForSTXXLSort {
				bool operator() (const packedBigPair& a, const packedBigPair& b) const {
					return a < b;
				}

				packedBigPair min_value() const { packedBigPair t; t.first = numeric_limits<uint64_t>::min(); t.second = numeric_limits<uint32_t>::min(); return t; }
				packedBigPair max_value() const { packedBigPair t; t.first = numeric_limits<uint64_t>::max(); t.second = numeric_limits<uint32_t>::max(); return t; }
			};

			stxxl::sort(tempPV.begin(), tempPV.end(), SCompareStructForSTXXLSort(), iMemorySizeAvail);

			Utilities::createFile(sDBPathOut);
			stxxlFile realOut(sDBPathOut, stxxl::file::RDWR);
			taxaOnly realIdx(&realOut, iSize);

			auto realIt = realIdx.begin();
			tOIt = tempPV.begin();
			while (tOIt != tempPV.end()) {
				*realIt = mIDsAsIdx.find(tOIt->second)->second;
				++realIt;
				++tOIt;
			}

			ofstream outSizeFile(sDBPathOut + "_info.txt");
			outSizeFile << tempPV.size();
			outSizeFile.close();
			ifstream oldFreqFile(indexFile + "_f.txt", std::ios::binary);
			ofstream newFreqFile(sDBPathOut + "_f.txt", std::ios::binary);

			newFreqFile << oldFreqFile.rdbuf();

			Trie T(static_cast<int8_t>(12), static_cast<int8_t>(iLowerK), iTrieDepth);
			T.SaveToStxxlVec(&tempPV, sDBPathOut);


			realIdx.export_files("_");
		}
		else {
			cout << "ERROR: No mode specified.\nPlease select a mode from the following list: generateCF, build, identify, shrink, update, delete, getFrequency, redundancy, trie." << endl;
		}
	} catch (const exception& e) {
		cerr << "ERROR: " << e.what() << endl;
		return 1;
	}

#if (_WIN32 || _WIN64) && _DEBUG
	system("PAUSE");
#endif
	return 0;
}