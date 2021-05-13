### TODO: Test for equality of output files (with cmp -s $old $new)
configfile: "snake_config.json"

rule all:
	input:
		json="work/results/example.json",
		json2="work/results/example_s.json",
		json3="work/results/example_u.json",
		freq="work/index/exampleIndex_u_f.txt",
		trie="work/index/exampleIndex_trie.txt",
		red="work/results/redundancy.txt",
		mergedIndex="work/index/index_merged",
		mergedIndex128="work/index/index_merged_128",
		pE="work/results/pairedEnd.json",
		jsonFromLarge="work/results/128.json",
		jsonFromProtein="work/results/protein.json",
		jsonFromTranslation="work/results/translate.json",
		idMult="work/results/mult.done"

rule generateCF:
	input:
		db="work/db/example.fasta",
		tax="taxonomy/",
		att="taxonomy/acc2tax/acc2Tax.txt"
	output:
		content="work/index/exampleIndex_content.txt"
	params:
		lvl="species",
		ram=config["ram"],
		callIdx="1"
	threads: config["threads"]
	shell:
		"{config[kASAExec]}kASA generateCF -c {output.content} -i {input.db} -f {input.att} -y {input.tax} -u {params.lvl} -n {threads} -m {params.ram} -v -x {params.callIdx}"

rule build:
	input:
		content=rules.generateCF.output.content,
		db="work/db/example.fasta"
	output:
		idx="work/index/exampleIndex",
		trieOut="work/index/exampleIndex_trie"
	params:
		ram=config["ram"],
		callIdx="2"
	threads: config["threads"]
	shell:
		"""
		{config[kASAExec]}kASA build -c {input.content} -d {output.idx} -i {input.db} -n {threads} -m {params.ram} -v -x {params.callIdx}
		cp work/index/exampleIndex_trie work/index/exampleIndex_trie_duplicate
		cp work/index/exampleIndex_trie.txt work/index/exampleIndex_trie_duplicate.txt
		"""

rule cloneContent:
	input:
		content=rules.generateCF.output.content
	output:
		contentOut="work/index/content_org.txt"
	shell:
		"cp {input.content} {output.contentOut}"

rule identify:
	input:
		idx=rules.build.output.idx,
		content=rules.cloneContent.output.contentOut,
		fq="work/input/example.fastq.gz"
	output:
		prof="work/results/example.csv",
		json="work/results/example.json"
	params:
		ram=config["ram"],
		callIdx="3"
	threads: config["threads"]
	shell:
		"{config[kASAExec]}kASA identify -c {input.content} -d {input.idx} -i {input.fq} -p {output.prof} -q {output.json} -n {threads} -m {params.ram} -v -x {params.callIdx}"

rule shrink:
	input:
		idx=rules.build.output.idx,
		content=rules.cloneContent.output.contentOut
	output:
		outIdx="work/index/exampleIndex_s"
	params:
		mode="2",
		ram=config["ram"],
		callIdx="4"
	threads: config["threads"]
	shell:
		"{config[kASAExec]}kASA shrink -c {input.content} -d {input.idx} -o {output.outIdx} -s {params.mode} -n {threads} -m {params.ram} -v -x {params.callIdx}"

rule identify_s:
	input:
		content=rules.cloneContent.output.contentOut,
		idx=rules.shrink.output.outIdx,
		fq="work/input/example.fastq.gz"
	output:
		prof="work/results/example_s.csv",
		json2="work/results/example_s.json"
	params:
		ram=config["ram"],
		callIdx="5"
	threads: config["threads"]
	shell:
		"{config[kASAExec]}kASA identify -c {input.content} -d {input.idx} -i {input.fq} -p {output.prof} -q {output.json2} -n {threads} -m {params.ram} -v -x {params.callIdx}"



rule update:
	input:
		dummyLink=rules.cloneContent.output.contentOut,
		db="work/db/16S_NCBI.fasta",
		tax="taxonomy/",
		att="taxonomy/acc2tax/acc2Tax.txt",
		content=rules.generateCF.output.content,
		idx=rules.build.output.idx
	output:
		index="work/index/exampleIndex_u"
	params:
		lvl="species",
		ram=config["ram"],
		callIdx="6"
	threads: config["threads"]
	shell:
		"""
		{config[kASAExec]}kASA update -d {input.idx} -i {input.db} -o {output.index} -f {input.att} -y {input.tax} -u {params.lvl} -n {threads} -m {params.ram} -v -x {params.callIdx}
		cp work/index/exampleIndex_u_f.txt work/index/exampleIndex_u_f_duplicate.txt
		"""

rule identify_u:
	input:
		idx=rules.update.output.index,
		fq="work/input/exampleInput.fasta"
	output:
		prof="work/results/example_u.csv",
		json3="work/results/example_u.json"
	params:
		ram=config["ram"],
		callIdx="7"
	threads: config["threads"]
	shell:
		"{config[kASAExec]}kASA identify -d {input.idx} -i {input.fq} -p {output.prof} -q {output.json3} -n {threads} -m {params.ram} -v -x {params.callIdx}"



rule reconstructFrequency:
	input:
		dummyLink=rules.identify_u.output.json3,
		index="work/index/exampleIndex_u"
	output:
		freqRec="work/index/exampleIndex_u_f.txt"
	params:
		ram=config["ram"],
		callIdx="8"
	threads: config["threads"]
	shell:
		"""
		{config[kASAExec]}kASA getFrequency -d {input.index} -n {threads} -m {params.ram} -v -x {params.callIdx}
		cmp -s work/index/exampleIndex_u_f_duplicate.txt work/index/exampleIndex_u_f.txt
		"""


rule reconstructTrie:
	input:
		linkDummy3="work/results/example.json",
		index=rules.build.output.idx,
		content=rules.cloneContent.output.contentOut
	output:
		trie="work/index/exampleIndex_trie.txt"
	params:
		ram=config["ram"],
		callIdx="9"
	threads: config["threads"]
	shell:
		"""
		{config[kASAExec]}kASA trie -c {input.content} -d {input.index} -n {threads} -m {params.ram} -v -x {params.callIdx}
		cmp -s work/index/exampleIndex_trie_duplicate work/index/exampleIndex_trie
		"""

rule checkRedundancy:
	input:
		dummyLink=rules.identify_u.output.json3,
		index=rules.update.output.index,
	output:
		red="work/results/redundancy.txt"
	params:
		ram = config["ram"],
		callIdx="10"
	threads: config["threads"]
	shell:
		"{config[kASAExec]}kASA redundancy -d {input.index} -n {threads} -m {params.ram} -v -x {params.callIdx} > {output.red}"

rule mergeIndices:
	input:
	output:
		mergedIndex="work/index/index_merged",
		mergedIndex128="work/index/index_merged_128"
	params:
		ram = config["ram"],
		callIdx="11"
	threads: config["threads"]
	shell:
		"""
		{config[kASAExec]}kASA build -i work/db/example.fasta -d work/index/index_1 -n {threads} -m {params.ram} -x {params.callIdx} -y taxonomy/ -f taxonomy/acc2tax/acc2Tax.txt -u "species"
		{config[kASAExec]}kASA build -i work/db/16S_NCBI.fasta -d work/index/index_2 -n {threads} -m {params.ram} -x {params.callIdx} -y taxonomy/ -f taxonomy/acc2tax/acc2Tax.txt -u "species"
		
		{config[kASAExec]}kASA merge -i work/index/index_1 -d work/index/index_2 -o work/index/index_merged -x {params.callIdx} -n {threads} -m {params.ram}
		
		{config[kASAExec]}kASA build -i work/db/example.fasta -d work/index/index_1_128 -n {threads} -m {params.ram} -x {params.callIdx} -y taxonomy/ -f taxonomy/acc2tax/acc2Tax.txt -u "species" --kH 25
		{config[kASAExec]}kASA build -i work/db/16S_NCBI.fasta -d work/index/index_2_128 -n {threads} -m {params.ram} -x {params.callIdx} -y taxonomy/ -f taxonomy/acc2tax/acc2Tax.txt -u "species" --kH 25
		
		{config[kASAExec]}kASA merge -i work/index/index_1_128 -d work/index/index_2_128 -o work/index/index_merged_128 -x {params.callIdx} -n {threads} -m {params.ram}
		"""

rule pairedEnd:
	input:
		dummyLink=rules.identify_u.output.json3,
		index=rules.update.output.index,
		freq=rules.reconstructFrequency.output.freqRec
	output:
		pE="work/results/pairedEnd.json"
	params:
		ram = config["ram"],
		callIdx="12"
	threads: config["threads"]
	shell:
		"{config[kASAExec]}kASA identify -d {input.index} -n {threads} -m {params.ram} -x {params.callIdx} -1 work/input/example.fastq.gz -2 work/input/example2.fastq.gz -q work/results/pairedEnd.json"

rule LargeK:
	input:
	output:
		jsonFromLarge="work/results/128.json"
	params:
		ram = config["ram"],
		callIdx="13"
	threads: config["threads"]
	shell:
		"""
		{config[kASAExec]}kASA build -i work/db/16S_NCBI.fasta -d work/index/index_128 -n {threads} -m {params.ram} -x {params.callIdx} -y taxonomy/ -f taxonomy/acc2tax/acc2Tax.txt -u "species" --kH 25
		{config[kASAExec]}kASA identify -d work/index/index_128 -n {threads} -m {params.ram} -x {params.callIdx} -i work/input/exampleInput.fasta -q work/results/128.json -k 25 7
		"""

rule Protein:
	input:
		content=rules.generateCF.output.content
	output:
		jsonFromProtein="work/results/protein.json"
	params:
		ram = config["ram"],
		callIdx="14"
	threads: config["threads"]
	shell:
		"""
		{config[kASAExec]}kASA build -c {input.content} -i work/db/ProtVulg_protein.fasta -d work/index/index_prot -n {threads} -m {params.ram} -x {params.callIdx} -z
		{config[kASAExec]}kASA identify -c {input.content} -d work/index/index_prot -n {threads} -m {params.ram} -x {params.callIdx} -z -i work/input/exampleProtein.fasta -q work/results/protein.json
		"""

rule ProteinWithAlphabet:
	input:
		content=rules.generateCF.output.content
	output:
		jsonFromTranslation="work/results/translate.json"
	params:
		ram = config["ram"],
		callIdx="15"
	threads: config["threads"]
	shell:
		"""
		{config[kASAExec]}kASA build -c {input.content} -i work/db/example.fasta -d work/index/index_translate -n {threads} -m {params.ram} -x {params.callIdx} -a work/table.prt 1
		{config[kASAExec]}kASA identify -c {input.content} -d work/index/index_translate -n {threads} -m {params.ram} -x {params.callIdx} -z -i work/input/exampleProtein.fasta -q work/results/translate.json
		"""

rule IdentifyMultiple:
	input:
		dummyLink=rules.identify_u.output.json3,
		idx=rules.update.output.index,
		derp=rules.reconstructFrequency.output.freqRec
	output:
		touch("work/results/mult.done")
	params:
		ram = config["ram"],
		callIdx="16"
	threads: config["threads"]
	shell:
		"""
		{config[kASAExec]}kASA identify_multiple -d {input.idx} -i work/input/ -q work/results/mult_ -n {threads} -m {params.ram} -x {params.callIdx}
		{config[kASAExec]}kASA identify_multiple -d {input.idx} -i work/input/ -q work/results/mult_r_ -n {threads} -m {params.ram} -x {params.callIdx} -r
		"""