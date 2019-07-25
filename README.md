# kASA

[![Current Version](https://img.shields.io/badge/version-1.0.7-green.svg)](https://github.com/SilvioWeging/kASA) [![GitHub issues](https://img.shields.io/github/issues/SilvioWeging/kASA.svg)](https://github.com/SilvioWeging/kASA/issues) ![GitHub All Releases](https://img.shields.io/github/downloads/SilvioWeging/kASA/total.svg)

This is the official repository of kASA - <u>k</u>-Mer <u>A</u>nalysis of <u>S</u>equences based on <u>A</u>mino acids, the preprint can be found [here](https://doi.org/10.1101/713966).

## Table of content
- [Things to know](#things-to-know-before-you-start)
- [Prerequisites](#prerequisites)
- [Setup](#setup)
	* [Linux](#linux)
	* [OS X](#OS X)
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
In order to do this, you need to set kASA up locally (no admin account needed!), create an index out of the genomic reference and then put in a file containing Reads.

Words like `<this>` are meant as placeholders to be filled with your specifics e.g. name, paths, ...

Folders and paths are recognized as such by letting a parameter end with a "/". 

## Prerequisites

Some scripts in the `/scripts` folder need Python 3.*, others are shell scripts. Most can be used just for convenience (see [Useful scripts](#useful-scripts)) but are not necessary for kASA.

You can use the system specific pre-compiled binaries in the `/bin` folder but I cannot guarantee that they will be universal enough to run on your system as well.

Note, that kASA is a console application so if you want to use these binaries, you must either use a terminal (Linux, OS X, Linux Subsystem for Windows) or PowerShell (Windows). A GUI may be implemented, depending on the amount of requests in the [poll](https://github.com/SilvioWeging/kASA/issues/1). If you're using the PowerShell, don't forget to add ".exe" at the end of each call to kASA: `.\<path to kASA>\kASA.exe`.

If you need to compile the code, you'll definitely need a C\+\+ compiler that supports C\+\+14 or if possible C\+\+17 (for `filesystem` and `execution`). This relates to Visual Studio 2017 and at least GCC version 6.1. On Linux and OS X, cmake is needed as well.

kASA depends on the [STXXL](https://stxxl.org/) and [Gzstream](https://www.cs.unc.edu/Research/compgeom/gzstream/) but contains all necessary files so you don't need to download those.

Last but not least: kASA provides an error (starts with "ERROR: ") and an output (starts with "OUT: ") stream. You can seperate them with 2> or 1>.

## Setup

### Linux

If you are running Linux in a Virtual Box and your host is Windows, you must do this first:
 
* Open Powershell in admin mode
* type `cd '<path to VirtualBox Installation>'`
* type `.\VBoxManage.exe setextradata <VM_NAME> VBoxInternal2/SharedFoldersEnableSymlinksCreate/<SHARE_NAME> 1`
* start virtual environment

If you are using a Linux distribution or the Linux Subsystem on Windows, you can proceed as follows:

Clone the repository.

Open a terminal and go to the path with `cd <installPath>/zlib/zlibBuild`.

Type `../configure` and after that `make`.

Now for kASA itself, type the following commands:

* `cd <installPath>`
* `mkdir <build folder name>`
* `cd <build folder name>`
* `cmake -DCMAKE_BUILD_TYPE=Release ..`
* `make`

### OS X

Since the native C\+\+ compiler on OS X (Clang) is not compatible with kASA(yet?), you'll need GCC (and cmake). Here are a few links, where you can get it:
* [Command Line Tools](http://osxdaily.com/2014/02/12/install-command-line-tools-mac-os-x/)
* [Homebrew](https://brew.sh/)
* [GCC](https://discussions.apple.com/thread/8336714)
* cmake: open a terminal and type `brew install cmake`

Afterwards (thanks for going through this trouble), please type this into your terminal:
```
	export CC=/usr/local/bin/gcc
	export CXX=/usr/local/bin/g++
```
and proceed as in the Linux part starting from "Clone ...".


### Windows

Clone the repository and open the file `slnForVS/kASA.sln` with Visual Studio 17.

Do [this](https://docs.microsoft.com/en-us/cpp/windows/how-to-use-the-windows-10-sdk-in-a-windows-desktop-application?view=vs-2017) to update to your Windows SDK-Version.

Switch to Release mode with x64.

Compile.

Change Parameters in Property Page &rarr; Debugging &rarr; Command Arguments.

You don't need to include "kASA" at the beginning like in the examples. Just start right with the mode e.g. `identify` when specifying parameters.

Run without debugging.

## TL;DR
```
<path to kASA>/kASA build -d <path and name of index file to be build> -i <fasta or folder with fastas> -m <amount of available GB> -n <number of parallel threads> -f <accToTaxFile(s)> -y <folder with nodes.dmp and names.dmp> -u <taxonomic level, e.g. species> <verbose>
e.g.: [weging@example ~] kASA/kASA build -d work/exampleIndex -i work/example.fasta -m 8 -n 4 -f taxonomy/acc2Tax/ -y taxonomy/ -u species -v

<path to kASA>/kASA shrink -d <path and name of index file> -s <1 or 2> (-g <percentage>)
e.g.: [weging@example ~] kASA/kASA shrink -d work/exampleIndex -s 2

<path to kASA>/kASA identify -d <path and name of small index file> -i <input file> -p <path and name of profile output> -q <path and name of read wise analysis> -m <amount of available GB> -n <number of parallel threads>
e.g.: [weging@example ~] kASA/kASA identify -d work/exampleIndex_s -i work/example.fastq.gz -p work/results/example.csv -q work/results/example.json -m 8 -n 4
```

## Modes and parameters

In this part, you can learn how to `build` your index, `identify` known stuff, `update` your index or `shrink` it to fit onto smaller platforms like external drives.

The first letter after a "-" is the short version of a parameter, "--" is the longer one. If no `<...>` follows a parameter, it's a boolean flag.

### Basic
Some parameters which are used by most modes:
##### Mandatory
* `-d (--database) <file>`: Actually path and name of the index but let's call it database since `-i` was already taken...
* `-c (--content) <file>`: Points to the content file. Can be defaulted when calling `build`, see [here](#generate-a-content-file).
* `-o (--outgoing) <file>`: Only important for `shrink`, `update` and `generateCF`. This file can be either the existing index or a new file name, depending on whether you want to keep the old file or not.

##### Optional
* `-t (--temp) <path>`: Path to temporary directory where files are stored that are either deleted automatically or can savely be deleted after kASA finishes. Defaults depend on your OS, [here](https://stxxl.org/tags/1.4.1/install_config.html) are some details.
* `-n (--threads) <number>`: Number of parallel threads. Note, that when compiling with C\+\+17 enabled on Windows, some routines use all available cores provided by the hardware (because of the implementation of parallel STL algorithms). Default: 1.
* `-m (--memory) <number>`: Amount of Gigabytes available to kASA. If you don't provide enough, a warning will be written and it attempts to use as little as possible but may crash. If you provide more than your system can handle, it will crash or thrash. If you write "inf" instead of a number, kASA assumes that you have no memory limit. Default: 5 GB.
* `-x (--callidx) <number>`: Number given to this call of kASA so that no problems with temporary files occur if multiple instances of kASA are running at the same time. Default: 0.
* `-v (--verbose)`: Prints out a little more information e.g. how much percent of your input was already read and analysed (if your input is not gzipped). Default: off.
* `-a (--alphabet) <file> <number>`: If you'd like to use a different translation alphabet formated in the NCBI compliant way, provide the file (gc.prt) and the id. Default: Hardcoded translation table.

### Generate a content file
##### Context
To fully use kASA, you first need a genomic database that you can create by either concatenating some fasta files (containing DNA) into one big file or putting them into a folder. From this data, a so called *content file* is created via the mode `generateCF`. This file contains a mapping of the accession numbers to the corresponding taxonomic IDs and names by using the NCBI taxonomy. Since it's only a text file, you can edit it if you want or just copy one from somewhere else. 

Furthermore, you'll need the NCBI taxonomy files `nodes.dmp` and `names.dmp` from [here](ftp://ftp.ncbi.nih.gov/pub/taxonomy/) contained in one of the `taxdump.*` files. 
The NCBI taxonomy offers multiple files for a mapping from accession number to taxid (see [here](ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/)) . 
If you know beforehand which one contains your mappings, just hand this to kASA. If not, please put them in one folder and hand the path to kASA instead. It's not necessary to uncompress them (the `.gz` at the end determines in which mode it'll be read).

Accession numbers in the fasta file(s) should be placed either right after the ">" e.g. ">CP023965.1 Proteus vulgaris" or inside the old format e.g. ">gi|89106884|ref|AC_000091.1| Escherichia coli str. K-12 substr. W3110". Anything else will get a dummy taxid.

If the content file contains entries with "EWAN_...", this stands for "Entries Without Accession Numbers". "unnamed" means they have a taxid but no name could be found (maybe due to deprecation).

This mode can be coupled with [Build](#build) by calling `build` and providing the same parameters described here but leaving out `-o` and `-c`. This creates a content file next to the index named `<index name>_content.txt` which then is considered the default. This eliminates the necessity of providing the `-c` parameter in almost every call.

##### Necessary parameters
* `-u (--level) <level>`: Taxonomic level at which you want to operate. Choose "lowest" if you want no linkage at a higher node in the taxonomic tree. All levels used in the NCBI taxonomy are available as well. To name a few: subspecies, species, genus, family, order, class, phylum, kingdom, superkingdom. 
* `-f (--acc2tax) <folder or file>`: As mentioned, either the folder containing the translation tables from accession number to taxid or a specific file.
* `-y (--taxonomy)` <folder>: This folder should contain the `nodes.dmp` and the `names.dmp` files.
* `-o (--outgoing) <file>`: Here, this parameter specifies where the content file should be written to.
##### Example call 
```
<path to kASA>/kASA generateCF -i <fastaFile(s)> -o <content file> -f <accToTaxFile(s)> -y <folder with nodes.dmp and names.dmp> -u <taxonomic level, e.g. species> (-v )
e.g.: [weging@example ~] kASA/kASA generateCF -i work/example.fasta -o work/content.txt -f taxonomy/acc2Tax/ -y taxonomy/ -u species -v
```

### Build
##### Context
This mode creates an index file, a frequency file(containing the amount of k-mers for each taxon) and a prefix trie out of the fasta file(s) you used in the previous mode.

This step can take much space and time, depending on the complexity and size of your database. Native support of translated sequences as a database will be added in a future update.

The content file from the previous mode is given to kASA via the `-c` parameter or can be created together with the index.

##### Optional paramameters
* `-c (--content) <file>`: Path and name of the content file either downloaded or created from genomic data.
* `-a (--alphabet) <file> <number>`: If you'd like to use a different translation alphabet formated in the NCBI compliant way, provide the file (gc.prt) and the id. Default: Hardcoded translation table.
##### Example call
```
<path to kASA>/kASA build -c <content file> -d <path and name of the index file> -i <folder or file> -t <temporary directory> -m <amount of RAM kASA can use> -n <number of threads>
e.g.: [weging@example ~] kASA/kASA build -c work/content.txt -d work/exampleIndex -i work/example.fasta -m 8 -t work/tmp/ -n 4

Create content file and index:
[weging@example ~] kASA/kASA build -d work/exampleIndex -i work/example.fasta -m 8 -t work/tmp/ -n 4 -f taxonomy/acc2Tax/ -y taxonomy/ -u species -v
```

### Identify
##### Context
This mode compares sequencing data with the built index. 

You can input fasta or fastq files, gzipped or not. It all depends on whether your file ends with `.gz` and if the first symbol is a `>` or `@`. If you want to put in multiple files, move them to a folder and place the path next to `-i`. The string given in `-p` and `-q` will then serve as a prefix concatenated with `"_<filename without path>.<csv or json>"`.

To input translated sequences, add the `-z` flag. If you've used a custom alphabet for conversion, just use the same here by copying the `-a <file> <number>` part of your `build` call.

Since kASA uses k-mers, a `k` can be given to influence accuracy. You can set these bounds by yourself with the `-k` parameter, the default lower bound is 7, the upper 12 (highest possible).
Smaller `k`'s than 6 only make sense if your data is very noisy or you're working on amino acid level.
If your read length is smaller than ![equation](http://www.sciweavers.org/tex2img.php?eq=%24k_%7Blower%7D%20%5Ccdot%203%24&bc=White&fc=Black&im=png&fs=12&ff=arev&edit=0) on DNA/RNA level, it will be padded. 

Another important thing here is the output. Or the output**s** if you want. kASA can give you two files, one contains the per-read information, which taxa were found (identification file, in json format) and the other a table of how much of each taxon was found (the profile, a csv file).
But because too much information isn't always nice, you can specify how much taxa shall be shown for each read and if the profile should be human readable. 


##### Necessary paramameters
* `-p (--profile) <file>`: Path and name of the profile that is put out.
* `-q (--rtt <file>)`: Path and name of the read ID to tax IDs output file. If not given, a profile-only version of kASA will be used which is much faster!
##### Optional paramameters
* `-z (--translated)`: Tell kASA, that the input consists of protein sequences. Note, that the index must've been
 converted via the same alphabet to amino acids.
* `-k <upper> <lower>`: Bounds for `k`, all `k`'s in between will be evaluated as well. If your intuition is more like `<lower> <upper>` then that's okay too. Default: 12 7.
* `--kH <upper>`: Set only the upper bound
* `--kL <lower>`: Set only the lower bound
* `-b (--beasts) <number>`: Number of hit taxa shown for each read. Default: 3.
* `-h (--human)`: Changes the output of the profile so that only necessary information is provided, see below.
##### Example call
```
<path to kASA>/kASA identify -c <content file> -d <path and name of index file> -i <input file or folder> -p <path and name of profile output> -q <path and name of read wise analysis> -m <amount of available GB> -t <path to temporary directory> -k <highest k> <lowest k> -n <number of parallel threads>
e.g.: [weging@example ~] kASA/kASA identify -c work/content.txt -d work/exampleIndex -i work/example.fastq.gz -p work/results/example.csv -q work/results/example.json -m 8 -t work/tmp/ -k 12 9 -n 4
```
#### Output
##### Normal:
###### Identification
```
[
	{
		"Read number": 0,
		"Specifier from input file": ">NZ_CP013542.1+NZ_JFYQ01000033.1",
		"Matched taxa": [
			{
				"tax ID": "396",
				"Name": "Rhizobium phaseoli",
				"k-mer Score": 33.32,
				"Relative Score": 1.03767e+00
			},
			{
				"tax ID": "1270",
				"Name": "Micrococcus luteus",
				"k-mer Score": 30.467,
				"Relative Score": 9.71243e-01
			}
		]
	}
]
```

###### Profile
|#tax ID|Name|Unique counts k=12|Unique counts k=11|...|Unique rel. freq. k=12|...|Non-unique counts k=12|...|Non-unique rel. freq. k=12|...)|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|9606,|Homo sapiens,|121252166,|111556464,|...|0.87,|...|2001658992,|...|0.79,|...|

##### Human readable:
###### Identification
|#Read number|Specifier from input file|Matched taxa|Names|Scores{relative,k-mer}|
|:---:|:---:|:---:|:---:|:---:|
|0 | @Dummy-line | 9606;147711 | Homo sapiens;Rhinovirus A | 1.332e+01,220.12;1.0221e+00,212.04 |

###### Profile
|tax ID|Name|Unique r.f. k=12|Non-unique r.f. k=12|
|:---:|:---:|:---:|:---:|
|9606,|Homo sapiens,|87%,|79%|

There are two relative frequencies because your index may be very ambigious (e.g. it only consists of E. Coli species) and thus has only few unique hits. 
To get a hint which one would be more relevant to you, check your index with a call to `redundancy` in the [Miscellaneous](#miscellaneous) section.

### Update
##### Context
Keeping the same index for years may not be that good an idea so kASA gives you the possibility to add genomic material to an existing index.
First, you need a fasta file or a folder with fasta files and a call to `update` with the `-o` parameter to specify where to put the new index if you don't want to overwrite the existing one.

Next, since the content file is updated as well, you'll need the same parameters as in [generateCF](#generate-a-content-file) (meaning `-u <...> -f <...> -y <...>`). If you've updated your content file manually then just add it via the `-c <...>` parameter.

If you want to delete entries from the index because they are not desired or deprecated in the NCBI taxonomy, add the `delnodes.dmp` file via the `-l` parameter.

It's not necessary to change the content file in this case although you should at some point to not clutter it too much...

If you've created the content file together with the index, this default content file will be used.

##### Necessary paramameters
* `-o (--outgoing) <file>`: Either the existing index or a new file name, depending on whether you want to keep the old file or not. Default: overwrite.
* `-l (--deleted) <file>`: delete taxa via the NCBI taxonomy file.
##### Example calls
```
<path to kASA>/kASA update -c <content file> -d <path and name of the index file> -o <path and name of the new index> -i <folder or file> -t <temporary directory> -m <amount of RAM> -f <accToTaxFile(s)> -y <folder with nodes.dmp and names.dmp> -u <taxonomic level, e.g. species>
e.g.: [weging@example ~] kASA/kASA identify -c work/content.txt -d work/exampleIndex -o work/updatedIndex -i work/newStuff.fasta -t work/tmp/ -m 8 -f taxonomy/acc2Tax/ -y taxonomy/ -u species

<path to kASA>/kASA delete -c <content file> -d <path and name of the index file> -o <path and name of the new index> -l <delnodes.dmp> -t <temporary directory> -m <amount of RAM>
e.g.: [weging@example ~] kASA/kASA identify -c work/content.txt -d work/exampleIndex -o work/updatedIndex -l taxonomy/delnodes.dmp -t work/tmp/ -m 8
``` 

### Shrink
##### Context 
Should updating your index not happen that often or you would like better performance and less space usage on your disk, shrinking it does the trick. kASA has two options:

1. The first way deletes a certain percentage of k-mers from every taxon. This may be lossy but impacts the accuracy not that much if your sequencing depth is high enough.
2. The second option is lossless as it assumes, that your content file is not larger than 65535 entries. This will reduce the size of the index by half but your index cannot be updated afterwars. Great for storing the index on an external drive to take with you.

The parameter `-o` also decides here, where to put your new index.
##### Necessary paramameters
* `-s (--strategy) <1 or 2>`: Shrink the index in the first or second way. This parameter may also be applied when building the index. Default is 2.
* `-g (--percentage) <integer>`: Deletes the given percentage of k-mers from every taxon. Default is 50.
* `-o (--outgoing) <file>`: Output path and name of your shrunken index file. Your other index cannot be overwritten with this. Default: takes your index file and appends a "_s".
##### Example call
```
<path to kASA>/kASA shrink -c <content file> -d <path and name of the index file> -o <path and name of the new index> -s <1 or 2> -g <percentage> -t <temporary directory>
e.g.: [weging@example ~] kASA/kASA shrink -c work/content.txt -d work/exampleIndex -o work/exampleIndex_s -s 2 -t work/tmp/
e.g.: [weging@example ~] kASA/kASA shrink -c work/content.txt -d work/exampleIndex -s 1 -g 25 -t work/tmp/
``` 

### Miscellaneous
1. If you've lost your frequency file "`<index name>_f.txt`", you can create one without building the index anew:
```
<path to kASA>/kASA getFrequency -c <content file> -d <path and name of the index file> -t <temporary directory> -m <amount of RAM> -n <number of threads>
``` 
2. If you've lost your trie file "`<index name>_trie`", your call's gonna look like:
```
<path to kASA>/kASA trie -d <path and name of the index file> -t <temporary directory>
``` 
If you've lost any of the other `.txt` files associated with your index: Bad luck! These contain the size of the data structure and cannot be restored, so you really need to re-create the corresponding file.

3. You can measure the redundancy of your index via:
```
<path to kASA>/kASA redundancy -c <content file> -d <path and name of the index file> -t <temporary directory>
```
This gives you a hint whether you should look at the unique relative frequencies and/or the non-unique relative frequencies. It's measured by counting how many tax IDs 99% of your k-mers have. If for example almost every k-mer has only one taxon, that's pretty unique.

4. The folder `example/` contains a minimum working example for most of the calls presented in this Readme with a dummy taxonomy, .fasta and .fastq.gz file.

## Useful scripts
- generateCF.py: The python 3 version of the `generateCF` mode described [here](#generate-a-content-file).
- jsonToFrequencies.py: Creates a profile based on the most prominent taxa per read. Useful, if your input could only be processed in many chunks (due to memory restrictions) which may influence the accuracy of the k-Mer profile negatively.

## Todos and upcoming
- A python script that creates Kraken-like output out of kASAs identification file
- Reworked building algorithm
- Join two built indices
- Profiles normalized to genome length, for now you could hack that with the frequency file
- New shrink mode deleting k-mers that are overrepresented
- Native support of the nr and other translated sequences
- Allow gzipped files as input for `build`
- Support of [Recentrifuge](https://github.com/khyox/recentrifuge)
- Support of [bioconda](https://bioconda.github.io/)/[Snakemake](https://snakemake.readthedocs.io/en/stable/)
- small collection of adapter sequences

## License

This project is licensed under the Boost License 1.0 - see the [LICENSE.md](LICENSE.md) file for details