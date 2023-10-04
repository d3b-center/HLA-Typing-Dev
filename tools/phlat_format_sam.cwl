cwlVersion: v1.2
class: CommandLineTool
id: phlat_format_sam
doc: "Format PHLAT bowtie2 output"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.threads)
  - class: DockerRequirement
    dockerPull: staphb/samtools:1.18

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      cut -f 3 $(inputs.phlat_sam.path) | grep -P "(A|B|C|DQA1|DQB1|DRB1)" > $(inputs.phlat_sam.nameroot)_HLA.map
  - position: 1
    shellQuote: false
    valueFrom: >-
      && awk  -F"\t" 'BEGIN {OFS="\t"} {if($3~/^A/){$4=$4-1+29910331;$8=$8-1+29910331}} {if($3~/^B/){$4=$4-1+31322260;$8=$8-1+31322260}}{if($3~/^C/){$4=$4-1+31236946;$8=$8-1+31236946}} {if($3~/^DQA1/){$4=$4-1+32605236;$8=$8-1+32605236}} {if($3~/^DQB1/){$4=$4-1+32628013;$8=$8-1+32628013}} {if($3~/^DRB1/){$4=$4-1+32546868;$8=$8-1+32546868}} {if($3~/^(A|B|C|D)/){$3="chr6"}} {$5="255"}{print $0}' $(inputs.phlat_sam.path)
  - position: 2
    shellQuote: false
    valueFrom: >-
      | cat $(inputs.bwt2_header.path) -
  - position: 3
    shellQuote: false
    valueFrom: >-
      | samtools view -b -o $(inputs.phlat_sam.nameroot).tmp.bam -@ $(inputs.threads)
  - position: 4
    shellQuote: false
    valueFrom: >-
      && samtools sort -n $(inputs.phlat_sam.nameroot).tmp.bam > $(inputs.phlat_sam.nameroot)_HLA.qnsorted.bam

inputs:
  phlat_sam: { type: File, doc: "Output sam file from bowtie2" }
  bwt2_header: { type: File, doc: "Desired header file to prepend to processed sorted bam"}
  threads: { type: 'int?', doc: "Num processing threads tp use", default: 8 }
  ram: { type: 'int?', doc: "Num GB memory to make available", default: 16 }

outputs:
  hla_map: { type: File, outputBinding: { glob: $(inputs.phlat_sam.nameroot)_HLA.map } }
  phlat_bam: { type: File, outputBinding: { glob: $(inputs.phlat_sam.nameroot).tmp.bam } }
  phlat_name_sort: { type: File, outputBinding: { glob: $(inputs.phlat_sam.nameroot)_HLA.qnsorted.bam } }
