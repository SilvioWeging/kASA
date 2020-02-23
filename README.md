# kASA

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/SilvioWeging/kASA) [![GitHub issues](https://img.shields.io/github/issues/SilvioWeging/kASA.svg)](https://github.com/SilvioWeging/kASA/issues) ![GitHub All Releases](https://img.shields.io/github/downloads/SilvioWeging/kASA/total.svg)

This is the official repository of kASA - <u>k</u>-Mer <u>A</u>nalysis of <u>S</u>equences based on <u>A</u>mino acids, the preprint can be found [here](https://doi.org/10.1101/713966).

## Table of content
- [Things to know](#things-to-know-before-you-start)
- [Prerequisites](#prerequisites)
- [Setup](#setup)
	* [Linux](#linux)
	* [macOS](#macOS)
	* [Windows](#windows)
- [TL;DR](#tl;dr)
- [Modes and paramameters](#modes-and-parameters)
	* [Basic](#basic)
	* [Content file](#generate-a-content-file)
	* [Build](#build)
	* [Identify](#identify)
		+ [Output](#output)
	* [Update](#update)
	* [Shrink](#shrink)
	* [Miscellaneous](#miscellaneous)
- [Useful scripts](#useful-scripts)
- [TODOS/Upcoming](#todos-and-upcoming)
- [License](#license)


## Things to know before you start

This tool is designed to read genomic sequences (also called Reads) and identify known parts by exactly matching k-mers to a reference database.
In order to do this, you need to set kASA up locally (no admin account needed!), create an index out of the genomic reference (for now, we will not provide standard indices) and then put in a file containing Reads.

Words like `<this>` are meant as placeholders to be filled with your specifics e.g. name, paths, ...

Folders and paths are recognized as such by letting a parameter end with a "/". 

## Prerequisites

Some scripts in the `/scripts` folder need Python 3.*, others are shell scripts. Most can be used just for convenience (see [Useful scripts](#useful-scripts)) but are not necessary for kASA.

You can use the system specific pre-compiled binaries in the `/bin` folder but I cannot guarantee that they will be universal enough to run on your system as well.

Note, that kASA is a console application so if you want to use these binaries, you must either use a terminal (Linux, macOS, Linux Subsystem for Windows) or PowerShell (Windows). A GUI may be implemented, depending on the amount of requests in the [poll](https://github.com/SilvioWeging/kASA/issues/1). If you're using the PowerShell, don't forget to add ".exe" at the end of each call to kASA: `.\<path to kASA>\kASA.exe`.

If you need to compile the code, you'll definitely need a C\+\+ compiler that supports C\+\+11 or if possible C\+\+17 (for `filesystem` and `execution`). We successfully tested it with Visual Studio 2017, GCC version 6.1, LLVM/Clang 9.0 and Apple Clang 9.0. [Here](https://en.cppreference.com/w/cpp/compiler_support) is a table showing if your version is sufficient. On Linux and macOS, cmake is needed as well.

kASA depends on the [STXXL](https://stxxl.org/) and [Gzstream](https://www.cs.unc.edu/Research/compgeom/gzstream/) but contains all necessary files so you don't need to download those.

Last but not least: kASA provides an error (starts with "ERROR: ") and an output (starts with "OUT: ") stream. You can seperate them with 2> or 1>.

## Setup

### Linux

Install cmake if you haven't already.

Open a terminal and navigate to the folder in which you would like to have kASA.

Clone the repository with `git clone https://github.com/SilvioWeging/kASA.git`.

First, build the zlib by going into the folder `zlib`.

Call `chmod +x configure` to give the file execution rights.

Create the folder `zlibBuild` and `cd` into it.

Type `../configure` and after that `make`.

Now for kASA itself, please type the following commands:

* `cd <installPath>/build` (or create/rename the folder)
* `cmake -DCMAKE_BUILD_TYPE=Release ..`
* `make`

### macOS

If you have Clang with LLVM installed (check with `clang --version`, is usually included in Xcode) do the same as above (Linux).

Should you prefer the GCC toolchain or don't have cmake installed, the following steps might be helpful:

* [Command Line Tools](http://osxdaily.com/2014/02/12/install-command-line-tools-mac-os-x/)
* [Homebrew](https://brew.sh/)
* [GCC](https://discussions.apple.com/thread/8336714)
* cmake: open a terminal and type `brew install cmake`

Afterwards, please type this into your terminal:
```
	export CC=/usr/local/bin/gcc
	export CXX=/usr/local/bin/g++
```
and proceed as in the Linux part starting from "Clone ...".


### Windows

Clone the repository and open the file `slnForVS/kASA.sln` with Visual Studio 17.

Do [this](https://docs.microsoft.com/en-us/cpp/windows/how-to-use-the-windows-10-sdk-in-a-windows-desktop-application?view=vs-2017) to update to your Windows SDK-Version.

Switch to Release mode with x64.

Build the project.

Change Parameters in Property Page &rarr; Debugging &rarr; Command Arguments.

You don't need to include "kASA" at the beginning like in the examples. Just start right with the mode e.g. `identify` when specifying parameters.

Run without debugging.

## TL;DR
```
build/kASA build -d <path and name of index file to be build> -i <fasta or folder with fastas> -m <amount of available GB> -n <number of parallel threads> -f <accToTaxFile(s)> -y <folder with nodes.dmp and names.dmp> -u <taxonomic level, e.g. species> <verbose>
e.g.: [weging@example:/kASA$] build/kASA build -d example/work/index/exampleIndex -i example/work/example.fasta -m 8 -n 2 -f example/taxonomy/acc2Tax/ -y example/taxonomy/ -u species -v

build/kASA identify -d <path and name of small index file> -i <input file> -p <path and name of profile output> -q <path and name of read wise analysis> -m <amount of available GB> -n <number of parallel threads>
e.g.: [weging@example:/kASA$] build/kASA identify -d example/work/index/exampleIndex -i example/work/example.fastq.gz -p example/work/results/example.csv -q example/work/results/example.json -m 5 -n 2
```

## Modes and parameters

In this part, you can learn how to `build` your index, `identify` known stuff, `update` your index or `shrink` it to fit onto smaller platforms like external drives.

The first letter after a "-" is the short version of a parameter, "--" is the longer one. If no `<...>` follows a parameter, it's a boolean flag.

### Basic
Some parameters which are used by most modes:
##### Mandatory
* `-d (--database) <file>`: Actually path and name of the index but let's call it database since `-i` was already taken...
* `-c (--content) <file>`: Points to the content file. Can be defaulted when calling `build`, see [here](#generate-a-content-file).
* `-o (--outgoing) <file>`: Only important for `shrink` and `update`. This file can be either the existing index or a new file name, depending on whether you want to keep the old file or not.

##### Optional
* `<mode> --help`: Shows all parameters for a given mode.
* `-t (--temp) <path>`: Path to temporary directory where files are stored that are either deleted automatically or can savely be deleted after kASA finishes. Defaults depend on your OS, [here](https://stxxl.org/tags/1.4.1/install_config.html) are some details. Typically, the path of the executable is used.
* `-n (--threads) <number>`: Number of parallel threads. Recommendation for different settings (due to I/O bottleneck): HDD: 1, SSD: 2-4, RAM disk: 2-?. Note, that when compiling with C\+\+17 enabled on Windows, some routines use all available cores provided by the hardware (because of the implementation of parallel STL algorithms). Default: 1.
* `-m (--memory) <number>`: Amount of RAM in Gigabytes available to kASA. If you don't provide enough, a warning will be written and it attempts to use as little as possible but may crash. If you provide more than your system can handle, it will crash or thrash. If you write "inf" instead of a number, kASA assumes that you have no memory limit. Default: 5 GB.
* `-x (--callidx) <number>`: Number given to this call of kASA so that no problems with temporary files occur if multiple instances of kASA are running at the same time. Default: 0.
* `-v (--verbose)`: Prints out a little more information e.g. how much percent of your input was already read and analysed (if your input is not gzipped). Default: off.
* `-a (--alphabet) <file> <number>`: If you'd like to use a different translation alphabet formated in the NCBI compliant way, provide the file (gc.prt) and the id (can be a string). Please use only letters in the range ['A',']'] from the ASCII table for your custom alphabet. Default: Hardcoded translation table.

### Generate a content file
##### Context
To fully use kASA, you first need a genomic database that you can create by either concatenating some fasta files (containing DNA) into one big file or putting them into a folder. From this data, a so called *content file* is created via the mode `generateCF`. This file contains a mapping of the accession numbers to the corresponding taxonomic IDs and names by using the NCBI taxonomy. Since it's only a text file, you can edit it if you want or just copy one from somewhere else. 

Furthermore, you'll need the NCBI taxonomy files `nodes.dmp` and `names.dmp` from [here](ftp://ftp.ncbi.nih.gov/pub/taxonomy/) contained in one of the `taxdump.*` files. 
The NCBI taxonomy offers multiple files for a mapping from accession number to taxid (see [here](ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/)) . 
If you know beforehand which one contains your mappings, just hand this to kASA. If not, please put them in one folder and hand the path to kASA instead. It's not necessary to uncompress them (the `.gz` at the end determines in which mode it'll be read).

Accession numbers in the fasta file(s) should be placed either right after the ">" e.g. ">CP023965.1 Proteus vulgaris" or inside the old format e.g. ">gi|89106884|ref|AC_000091.1| Escherichia coli str. K-12 substr. W3110". Anything else will get a dummy taxid.

If the content file contains entries with "EWAN_...", this stands for "Entries Without Accession Numbers". "unnamed" means they have a taxid but no name could be found (maybe due to deprecation).

This mode can be coupled with [Build](#build) by calling `build` and providing the same parameters described here but leaving out `-c`. This creates a content file next to the index named `<index name>_content.txt` which then is considered the default. This eliminates the necessity of providing the `-c` parameter in almost every call.

Note, that this step is optional if you provide your own content file or another index with the same content file shall be created (this means creating subsets of a full index is possible with the same content file). The accepted format per line is as follows:
```
<Name>	<taxid of specified level e.g. species>	<taxids on lowest level>	<accession numbers>
```
which for example could look like this:
```
Proteus vulgaris	585	585	CP023965.1;NZ_NBUT01000031.1
Hyphomicrobium denitrificans	53399	582899;670307	NC_014313.1;NC_021172.1
```
The taxids are required to be integers and no header line is necessary. This file can be given to `build` via the `-c` parameter.

##### Necessary parameters
* `-i (--input) <file/folder>`: Fasta file(s). Can be gzipped (but must end with `.gz`). If you want to process multiple files at once, put them inside a folder and let the path end with `/`. No default.
* `-u (--level) <level>`: Taxonomic level at which you want to operate. All levels used in the NCBI taxonomy are available as well. To name a few: subspecies, species, genus, family, order, class, phylum, kingdom, superkingdom. Choose "lowest" if you want no linkage at a higher node in the taxonomic tree, this corresponds to other tools' "sequence" level. That means that no real taxid will be given and the name will be the line from the fasta containing the accession number. Default: species.
* `-f (--acc2tax) <folder or file>`: As mentioned, either the folder containing the translation tables from accession number to taxid or a specific file. Can be gzipped.
* `-y (--taxonomy)` <folder>: This folder should contain the `nodes.dmp` and the `names.dmp` files.
* `-c (--content) <file>`: Here, this parameter specifies where the content file should be written to.
##### Example call 
```
<path to kASA>/kASA generateCF -i <fastaFile(s)> -c <content file> -f <accToTaxFile(s)> -y <folder with nodes.dmp and names.dmp> -u <taxonomic level, e.g. species> (-v )
e.g.: [weging@example:/kASA$] build/kASA generateCF -i work/example.fasta -c work/content.txt -f taxonomy/acc2Tax/ -y taxonomy/ -u species -v
```

### Build
##### Context
This mode creates an index file, a frequency file (containing the amount of k-mers for each taxon) and a prefix trie out of the fasta file(s) you used in the previous mode.

This step can take much space and time, depending on the complexity and size of your database. Native support of translated sequences as a database will be added in a future update.

The content file from the previous mode is given to kASA via the `-c` parameter or can be created together with the index.

##### Necessary parameters
* `-i (--input) <file/folder>`: Fasta file(s). Can be gzipped (but must end with `.gz`). If you want to process multiple files at once, put them inside a folder and let the path end with `/`. No default.
* `-d (--database) <file>`: Actually path and name of the index but let's call it database since `-i` was already taken...
##### Optional paramameters
* `-c (--content) <file>`: Path and name of the content file either downloaded or created from genomic data.
* `-a (--alphabet) <file> <number>`: If you'd like to use a different translation alphabet formated in the NCBI compliant way, provide the file (gc.prt) and the id (can be a string). Please use only letters in the range ['A',']'] from the ASCII table for your custom alphabet. Default: Hardcoded translation table.
* `-z (--translated)`: Tell kASA, that the input consists of protein sequences. Currently in BETA.
##### Example call
```
<path to kASA>/kASA build -c <content file> -d <path and name of the index file> -i <folder or file> -t <temporary directory> -m <amount of RAM kASA can use> -n <number of threads>
e.g.: [weging@example:/kASA$] build/kASA build -c example/work/content.txt -d  example/work/index/exampleIndex -i example/work/example.fasta -m 8 -t example/work/tmp/ -n 2

Create content file and index:
[weging@example:/kASA$] build/kASA build -d  example/work/index/exampleIndex -i example/work/example.fasta -m 8 -t example/work/tmp/ -n 2 -f taxonomy/acc2Tax/ -y taxonomy/ -u species -v
```

### Identify
##### Context
This mode compares sequencing data with the built index. 

You can input fasta or fastq files, gzipped or not. It all depends on whether your file ends with `.gz` and if the first symbol is a `>` or `@`. If you want to put in multiple files, move them to a folder and place the path next to `-i`. The string given in `-p` and `-q` will then serve as a prefix concatenated with `"_<filename without path>.<csv or json>"`.

To input translated sequences, add the `-z` flag. If you've used a custom alphabet for conversion, just use the same here by copying the `-a <file> <number>` part of your `build` call.

Since kASA uses k-mers, a `k` can be given to influence accuracy. You can set these bounds by yourself with the `-k` parameter, the default lower bound is 7, the upper 12 (highest possible).
Smaller `k`'s than 6 only make sense if your data is very noisy or you're working on amino acid level.
If your read length is smaller than ![equation](http://www.sciweavers.org/tex2img.php?eq=%24k_%7Blower%7D%20%5Ccdot%203%24&bc=White&fc=Black&im=png&fs=12&ff=arev&edit=0) on DNA/RNA level, it will be padded. 

If you want to optimise precision over sensitivity, you could use `k 12 12` and/or filter out low scoring reads (e.g. by ignoring everything below 0.5 (Relative Score)).

Another important thing here is the output. Or the output**s** if you want. kASA can give you two files, one contains the per-read information, which taxa were found (identification file, in json format) and the other a table of how much of each taxon was found (the profile, a csv file).
But because too much information isn't always nice, you can specify how much taxa shall be shown for each read (e.g. `-b 5` shows best 5 hits) or if the outputs should be human readable (see below). 

The per read error score ranges from -1 to 1. A 1 means that the best score deviates as far as possible from the optimal score, 0 means a perfect match and -1 means that the reverse complement also fits perfectly. In human readable format, only the error of the best score is printed.

Note, that if you input a folder, file names are appended to your string given via `-p` or `-q`. If for example a folder contains two files named `example1.fq` and `example2.fasta` with `-p example/work/results/out_` as a parameter, then kASA will generate two output files named `out_example1.fq.csv` and `out_example2.fasta.csv`.

If a read cannot be identified, the array "Top hits" in json format is empty, and in human readable format "-" is printed in every column instead of taxa, names and scores. 
The "Top hits" array will only contain multiple entries, if their k-mer Score is identical. Otherwise it contains the entry with the highest relative Score and all other hits are saved into the "Further hits" array.

Should you provide more than 9 GB of RAM and a lower `k` of at least 7 (default), a much faster lookup table instead of a prefix trie is used.

The first line of the profile is always "not identified" followed by zeroes for the unique and non-unique frequencies but with values for the overall frequencies describing the fracture of the k-mers from the input, which could not be identified.

##### Necessary paramameters
* `-i (--input) <file/folder>`: Fastq or fasta file(s). Can be gzipped (but must end with `.gz`). If you want to process multiple files at once, put them inside a folder and let the path end with `/`. No default.
* `-p (--profile) <file>`: Path and name of the profile that is put out.
* `-q (--rtt) <file>`: Path and name of the read ID to tax IDs output file. If not given, a profile-only version of kASA will be used which is much faster!
##### Optional paramameters
* `-r (--ram)`: Loads the index into primary memory. If you don't provide enough RAM for this, it will fall back to using secondary memory. Default: false.
* `-z (--translated)`: Tell kASA, that the input consists of protein sequences. Note, that the index must've been
 converted via the same alphabet to amino acids.
* `-a (--alphabet) <file> <number>`: If you'd like to use a different translation alphabet formated in the NCBI compliant way, provide the file (gc.prt) and the id (can be a string). Please use only letters in the range ['A',']'] from the ASCII table for your custom alphabet. Default: Hardcoded translation table.
* `-k <upper> <lower>`: Bounds for `k`, all `k`'s in between will be evaluated as well. If your intuition is more like `<lower> <upper>` then that's okay too. Default: 12 7.
* `--kH <upper>`: Set only the upper bound
* `--kL <lower>`: Set only the lower bound
* `-b (--beasts) <number>`: Number of hit taxa shown for each read. Default: 3.
* `-h (--human)`: Changes the output of the profile so that only necessary information is provided, see below.
* `-e (--unique)`: Ignores duplicates of `k`-mers in every read. This helps removing bias from repeats but messes with error scores, so consider it BETA.
##### Example call
```
<path to kASA>/kASA identify -c <content file> -d <path and name of index file> -i <input file or folder> -p <path and name of profile output> -q <path and name of read wise analysis> -m <amount of available GB> -t <path to temporary directory> -k <highest k> <lowest k> -n <number of parallel threads>
e.g.: [weging@example:/kASA$] build/kASA identify -c example/work/content.txt -d  example/work/index/exampleIndex -i example/work/example.fastq.gz -p example/work/results/example.csv -q example/work/results/example.json -m 8 -t example/work/tmp/ -k 12 9 -n 2
```
#### Output
##### Normal:
###### Identification
```
[
{
	"Read number": 0,
	"Specifier from input file": "NZ_CP013542.1+NZ_JFYQ01000033.1",
	"Top hits": [
	{
		"tax ID": "396",
		"Name": "Rhizobium phaseoli",
		"k-mer Score": 33.32,
		"Relative Score": 1.03767e+00,
		"Error": 0.67
	}
	],
	"Further hits": [
	{
		"tax ID": "1270",
		"Name": "Micrococcus luteus",
		"k-mer Score": 31.42,
		"Relative Score": 9.71243e-01
		"Error": 0.69
	},
	{
		"tax ID": "2128",
		"Name": "Mycoplasma flocculare",
		"k-mer Score": 0.0833333,
		"Relative Score": 4.561024e-03,
		"Error": 0.99729
	}
	]
}
]
```

###### Profile
|#tax ID|Name|Unique counts k=12|Unique counts k=11|...|Unique rel. freq. k=12|...|Non-unique counts k=12|...|Non-unique rel. freq. k=12|...|Overall rel. freq. k=12| ... |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|9606,|Homo sapiens,|121252166,|111556464,|...|0.87,|...|2001658992,|...|0.79,|...|0.65,|...|

##### Human readable:
###### Identification (tab-separated)
|#Read number|Specifier from input file|Matched taxa|Names|Scores{relative,k-mer}|Error|
|:---:|:---:|:---:|:---:|:---:|:---:|
|0 | Dummy-line | 9606;147711 | Homo sapiens;Rhinovirus A | 1.332e+01,220.12;1.0221e+00,212.04 | -0.001 |

###### Profile
|#tax ID|Name|Unique rel. freq. in %|Non-unique rel. freq. in %|Overall rel. freq. in %|
|:---:|:---:|:---:|:---:|:---:|
|9606,|Homo sapiens,|87,|79,|65|

There are two relative frequencies because your index may be very ambigious (e.g. it only consists of E. Coli species) and thus has only few unique hits. 
To get a hint which one would be more relevant to you, check your index with a call to `redundancy` in the [Miscellaneous](#miscellaneous) section.

Relative frequencies in the human readable profile are given for the largest k. The "Overall relative frequency" is calculated by dividing the non-unique counts by the total number of k-mers from the input.
This is also printed in the verbose mode like: "OUT: Number of k-mers in input: ... of which ... % were identified." for the largest k.

### Update
##### Context
Keeping the same index for years may not be that good an idea so kASA gives you the possibility to add genomic material to an existing index.
First, you need a fasta file or a folder with fasta files and a call to `update` with the `-o` parameter to specify where to put the new index if you don't want to overwrite the existing one.

Next, since the content file is updated as well, you'll need the same parameters as in [generateCF](#generate-a-content-file) (meaning `-u <...> -f <...> -y <...>`). If you've updated your content file manually then just add it via the `-c <...>` parameter.

If you want to delete entries from the index because they are not desired or deprecated in the NCBI taxonomy, add the `delnodes.dmp` file via the `-l` parameter.
It's not necessary to change the content file in this case although you should at some point to not clutter it too much...

If you've created the content file together with the index, this default content file will be used.

##### Necessary paramameters
* `-i (--input) <file/folder>`: Fasta file(s). Can be gzipped (but must end with `.gz`). If you want to process multiple files at once, put them inside a folder and let the path end with `/`. No default.
* `-o (--outgoing) <file>`: Either the existing index or a new file name, depending on whether you want to keep the old file or not. Default: overwrite.
* `-l (--deleted) <file>`: delete taxa via the NCBI taxonomy file.
##### Optional paramameters
* `-z (--translated)`: Tell kASA, that the input consists of protein sequences. Note, that the index must've been
 converted via the same alphabet to amino acids. Currently in BETA.
* `-a (--alphabet) <file> <number>`: If you'd like to use a different translation alphabet formated in the NCBI compliant way, provide the file (gc.prt) and the id (can be a string). Please use only letters in the range ['A',']'] from the ASCII table for your custom alphabet. Default: Hardcoded translation table.
##### Example calls
```
<path to kASA>/kASA update -d <path and name of the index file> -o <path and name of the new index> -i <folder or file> -t <temporary directory> -m <amount of RAM> -f <accToTaxFile(s)> -y <folder with nodes.dmp and names.dmp> -u <taxonomic level, e.g. species>
e.g.: [weging@example:/kASA$] build/kASA update -c example/work/content.txt -d  example/work/index/exampleIndex -o example/work/updatedIndex -i example/work/newStuff.fasta -t example/work/tmp/ -m 8 -f taxonomy/acc2Tax/ -y taxonomy/ -u species

<path to kASA>/kASA delete -c <content file> -d <path and name of the index file> -o <path and name of the new index> -l <delnodes.dmp> -t <temporary directory> -m <amount of RAM>
e.g.: [weging@example:/kASA$] build/kASA delete -c example/work/content.txt -d  example/work/index/exampleIndex -o example/work/updatedIndex -l taxonomy/delnodes.dmp -t example/work/tmp/ -m 8
``` 

### Shrink
##### Context 
Should updating your index not happen that often or you would like better performance and less space usage on your disk, shrinking it does the trick. kASA has multiple options:

1. The first way deletes a certain percentage of k-mers from every taxon. This may be lossy but impacts the accuracy not that much if your sequencing depth is high enough.
2. The second option is lossless but it assumes, that your content file is not larger than 65535 entries and that you don't need k's smaller than 7. This will reduce the size of the index by half but your index cannot be updated afterwards. Great for storing the index on an external drive.
3. This lossy option determines the (normalized binary) entropy of every k-mer and throws away anything not containing enough information. For example: AAABBBAAABBB would be thrown away but ABCDEFGAAABC wouldn't.

The parameter `-o` also decides here, where to put your new index.
##### Necessary paramameters
* `-s (--strategy) <1, 2 or 3>`: Shrink the index in the first or second way. Default is 2.
* `-g (--percentage) <integer>`: Deletes the given percentage of k-mers from every taxon. This parameter may also be applied when building the index (for example: -g 50 skips every second k-mer).
* `-o (--outgoing) <file>`: Output path and name of your shrunken index file. Your other index cannot be overwritten with this. Default: takes your index file and appends a "_s".
##### Example call
```
<path to kASA>/kASA shrink -c <content file> -d <path and name of the index file> -o <path and name of the new index> -s <1 or 2> -g <percentage> -t <temporary directory>
e.g.: [weging@example:/kASA$] build/kASA shrink -c example/work/content.txt -d  example/work/index/exampleIndex -o  example/work/index/exampleIndex_s -s 2 -t example/work/tmp/
e.g.: [weging@example:/kASA$] build/kASA shrink -c example/work/content.txt -d  example/work/index/exampleIndex -s 1 -g 25 -t example/work/tmp/
e.g.: [weging@example:/kASA$] build/kASA build -c example/work/content.txt -d  example/work/index/exampleIndex -g 50 -i example/work/example.fasta -m 8 -t example/work/tmp/ -n 2
``` 

### Miscellaneous
1. If you've lost your frequency file "`<index name>_f.txt`" and have not shrunken the index via mode 2, you can create one without building the index anew:
```
<path to kASA>/kASA getFrequency -c <content file> -d <path and name of the index file> -t <temporary directory> -m <amount of RAM> -n <number of threads>
``` 

2. If you've lost your trie file "`<index name>_trie`" and have not shrunken the index via mode 2, your call's gonna look like:
```
<path to kASA>/kASA trie -d <path and name of the index file> -t <temporary directory>
``` 

3. If you've lost the `_info.txt` file associated with your index: Get the size in bytes, divide by 12 (not shrunken) or 6 (shrunken via method 2) and create a `<index name>_info.txt` file with the result as first line in it.
For the shrunken index, add a 3 as second line in the file (which indicates, that it has been shrunken).

4. You can measure the redundancy of your (non-halved) index via:
```
<path to kASA>/kASA redundancy -c <content file> -d <path and name of the index file> -t <temporary directory>
```
This gives you a hint whether you should look at the unique relative frequencies and/or the non-unique relative frequencies. It's measured by counting how many tax IDs 99% of your k-mers have. If for example almost every k-mer has only one associated taxon, that's pretty unique.

5. The folder `example/` contains a minimum working example for most of the calls presented in this Readme with a dummy taxonomy, .fasta and .fastq.gz file. If you have [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) installed, you can use the Snakemake file inside the folder to run all possible modes. Go inside the folder with `cd` and type `snakemake --config path=../build/` to point to the path where executable of kASA lies. If anything goes awry, please tell me.

## Useful scripts
- jsonToFrequencies.py: Creates a profile based on the most prominent taxa per read. Usage: `-i <kASA output> -o <result> (-t <threshold for rel. score>)`. Consumes a lot of memory because the json file is loaded into memory.
- hrToFrequencies.py: Same as above but for human readable output.
- csvToCAMI.py: Converts a profiling output into the CAMI profile format. Needs the NCBI 'nodes.dmp' and 'names.dmp' file for upwards traversing of the taxonomic tree. Usage: `-i <kASA output> -o <result> -n <nodes.dmp> -m <names.dmp> -k <k value> -u <u: unique, o: overall, n: non-unique> (-t <threshold>)`.
- camiToKrona.py: Converts the CAMI profile to a file format readable by [Krona](https://github.com/marbl/Krona/wiki). Usage: `-i <cami file> -o <krona file>`.
- jsonToCAMIBin.py: Converts the json output file into the CAMI binning format. Usage: `-i <json file> -o <cami file>`.


## Todos and upcoming
- A python script that creates Kraken-like output out of kASAs identification file
- ~~Reworked building algorithm~~
- Join two built indices
- Profiles normalized to genome length, for now you could hack that with the frequency file
- ~~New shrink mode deleting k-mers that are overrepresented~~
- ~~Native support of the nr and other translated sequences~~
- ~~Allow gzipped files as input for `build`~~
- ~~RAM mode~~
- ~~Omission of the prefix trie if enough RAM (ca 12GB) is available~~
- ~~Support of Clang (macOS)~~
- ~~Snakemake pipeline for quality control~~
- Support of [Recentrifuge](https://github.com/khyox/recentrifuge)
- Support of [bioconda](https://bioconda.github.io/)/[Snakemake](https://snakemake.readthedocs.io/en/stable/)
- small collection of adapter sequences
- consideration of paired-end information
- bzip2 support

## License

This project is licensed under the Boost License 1.0 - see the [LICENSE](https://github.com/SilvioWeging/kASA/blob/master/LICENSE.txt) file for details
