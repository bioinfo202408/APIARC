#!/usr/bin/perl
use strict;
use warnings;
use YAML 'LoadFile';
use Getopt::Long;

my ($config_file, $outdir, $picture_dir);
GetOptions(
    'config=s'  => \$config_file,
    'outdir=s'  => \$outdir,
    'picdir=s'  => \$picture_dir,
);

my %genome_to_annotation = (
    'mm' => 'reference/annotation/gencode.vM25.annotation.bed',
    'hg' => 'reference/annotation/gencode.v44.annotation.bed',
);

my $config = LoadFile($config_file);

system("mkdir -p $picture_dir");
foreach my $sample (keys %{$config}) {
    my $group = $config->{$sample};
    next unless ref($group) eq 'HASH';
    my ($treat_bam, $control_bam, $bw1, $bw2, $genome, $qval) = @$group{qw/treat_bam control_bam bw1 bw2 genome qval/};
    $qval //= 0.05;
    next unless $treat_bam && $control_bam && $bw1 && $bw2 && $genome;
    next unless exists $genome_to_annotation{$genome};
    my $annotation_file = $genome_to_annotation{$genome};
    my $sample_dir = "$outdir/$sample";
    system("mkdir -p $sample_dir");
    my $bam_control_dir = "$sample_dir/bam_control";
    system("mkdir -p $bam_control_dir");
    my $spp_pdf = "$picture_dir/${sample}_bam_control.q${qval}.pdf";
    my $spp_out = "$bam_control_dir/${sample}_cross.txt";
    my $spp_cmd = "run_spp.R -rf -c=$treat_bam -i=$control_bam -p=8 -odir=$bam_control_dir -savp=$spp_pdf -out=$spp_out";
    system($spp_cmd);
    system("macs2 callpeak -t $treat_bam -c $control_bam -g $genome -n $sample --keep-dup all -q $qval --outdir $sample_dir");
    my $narrow = "$sample_dir/${sample}_peaks.narrowPeak";
    my $bed = "$sample_dir/${sample}_peaks.bed";
    my $tab_bed = "$sample_dir/${sample}_peaks_tab.bed";
    system("awk '{print \$1, \$2, \$3, \$4, \$5}' $narrow > $bed");
    system("sed 's/ \\+/\t/g' $bed > $tab_bed");
}
my %sample_groups;
foreach my $sample (keys %{$config}) {
    my $group = $config->{$sample};
    next unless ref($group) eq 'HASH';
    my ($sample_type) = $sample =~ /^(.*)_rep\d+$/;
    $sample_type ||= $sample;
    push @{$sample_groups{$sample_type}}, $group->{bw1} if $group->{bw1};
    push @{$sample_groups{$sample_type}}, $group->{bw2} if $group->{bw2};
}
my @group_matrix_pairs;
foreach my $sample_type (keys %sample_groups) {
    my @sample_names = grep { /^$sample_type\_rep\d+$/ && ref($config->{$_}) eq 'HASH' } keys %$config;
    @sample_names = sort @sample_names;
    next unless @sample_names;
    my @bw_files = map { $config->{$_}->{bw1} } @sample_names;
    my $input_bw = $config->{$sample_names[0]}->{bw2};
    push @bw_files, $input_bw if $input_bw;
    my @labels = (@sample_names, 'Input');
    my $bw_files_str = join(" ", @bw_files);
    my $labels_str = join(" ", @labels);
    my $output_dir = "$outdir/$sample_type";
    system("mkdir -p $output_dir");
    my $rep1_sample = "${sample_type}_rep1";
    my $genome = $config->{$rep1_sample} ? $config->{$rep1_sample}->{genome} : undef;
    next unless $genome && exists $genome_to_annotation{$genome};
    my $annotation_file = $genome_to_annotation{$genome};
    my $matrix_output = "$output_dir/${sample_type}_matrix.gz";
    my $matrix_values = "$output_dir/${sample_type}_matrix_values.tab";
    my $compute_matrix_cmd = "computeMatrix reference-point -R $annotation_file -S $bw_files_str -b 1000 -a 1000 --binSize 10 -p 20 -o $matrix_output --outFileNameMatrix $matrix_values --samplesLabel $labels_str";
    system($compute_matrix_cmd);
    my $sample_names_str = join(",", @sample_names);
    my $input_name = 'Input';
    push @group_matrix_pairs, join("\t", $sample_type, $matrix_values, $sample_names_str, $input_name);
}
my $txt_file = "$outdir/group_matrix.txt";
open my $txt_fh, ">", $txt_file or die "Cannot open $txt_file: $!";
foreach my $pair (@group_matrix_pairs) {
    print $txt_fh "$pair\n";
}
close $txt_fh;

# nohup perl peakcalling_script.pl --config /home/yxiaobo/RCProj/ChIPseq/database/TF-binding/GSE85632_text/bambw_list.yaml --outdir /home/yxiaobo/RCProj/ChIPseq/database/TF-binding/GSE85632_text/Macs2 --picdir /home/yxiaobo/RCProj/picture/ChIPseq/TF-binding/GSE85632_text  &
