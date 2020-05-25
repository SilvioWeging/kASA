### TODO: Test for equality of output files (with cmp -s $old $new), all modes, protein input, protein database, different alphabet, 
configfile: "snake_config.json"

rule all:
	input:
		json="work/results/example.json",
		json2="work/results/example_s.json",
		json3="work/results/example_u.json",
		freq="work/index/exampleIndex_u_f.txt",
		trie="work/index/exampleIndex_trie.txt",
		red="work/results/redundancy.txt"

rule generateCF:
	input:
		db="work/example.fasta",
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
		db="work/example.fasta"
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
		fq="work/example.fastq.gz"
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
		fq="work/example.fastq.gz"
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
		db="work/16S_NCBI.fasta",
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
		{config[kASAExec]}kASA update -c {input.content} -d {input.idx} -i {input.db} -o {output.index} -f {input.att} -y {input.tax} -u {params.lvl} -n {threads} -m {params.ram} -v -x {params.callIdx}
		cp work/index/exampleIndex_u_f.txt work/index/exampleIndex_u_f_duplicate.txt
		"""

rule identify_u:
	input:
		content="work/index/exampleIndex_content.txt",
		idx=rules.update.output.index,
		fq="work/exampleInput.fasta"
	output:
		prof="work/results/example_u.csv",
		json3="work/results/example_u.json"
	params:
		ram=config["ram"],
		callIdx="7"
	threads: config["threads"]
	shell:
		"{config[kASAExec]}kASA identify -c {input.content} -d {input.idx} -i {input.fq} -p {output.prof} -q {output.json3} -n {threads} -m {params.ram} -v -x {params.callIdx}"



rule reconstructFrequency:
	input:
		index="work/index/exampleIndex_u",
		content="work/index/exampleIndex_content.txt"
	output:
		freqRec="work/index/exampleIndex_u_f.txt"
	params:
		ram=config["ram"],
		callIdx="8"
	threads: config["threads"]
	shell:
		"""
		{config[kASAExec]}kASA getFrequency -c {input.content} -d {input.index} -n {threads} -m {params.ram} -v -x {params.callIdx}
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
		index=rules.update.output.index,
		content=rules.generateCF.output.content
	output:
		red="work/results/redundancy.txt"
	params:
		ram = config["ram"],
		callIdx="10"
	threads: config["threads"]
	shell:
		"{config[kASAExec]}kASA redundancy -c {input.content} -d {input.index} -n {threads} -m {params.ram} -v -x {params.callIdx} > {output.red}"