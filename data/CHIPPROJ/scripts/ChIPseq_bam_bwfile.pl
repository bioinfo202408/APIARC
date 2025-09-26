#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use File::Basename;
use POSIX qw(waitpid);

my ($inputdir, $outputdir, $indexdir, $picarddir, $threads, $help, $single_end);
GetOptions(
    "inputdir|i=s"   => \$inputdir,
    "outputdir|o=s"  => \$outputdir,
    "indexdir|x=s"   => \$indexdir,
    "picarddir|p=s"  => \$picarddir,
    "threads|t=i"    => \$threads,
    "single_end|s!"  => \$single_end,
    "help!"          => \$help,
);

die "Usage: perl readsmapping.pl --inputdir <dir> --outputdir <dir> --indexdir <dir> --picarddir <dir> --threads <num> [--single_end]\n"
    if $help || !$inputdir || !$outputdir || !$indexdir || !$picarddir || !$threads;

make_path($outputdir) unless -d $outputdir;
my $bowtie_log_all = "$outputdir/bowtie.log";

my @supported_formats = ("fastq", "fq", "fastq.gz", "fq.gz");
my @fq1_files;
foreach my $format (@supported_formats) {
    push @fq1_files, glob("$inputdir/*_1.$format");
}

my %sample_info;
foreach my $fq1 (@fq1_files) {
    if ($fq1 =~ /(.+)_1\.(fastq|fq|fastq\.gz|fq\.gz)$/) {
        my $base = $1;
        my $extension = $2;
        my $sample_id = basename($base);
        my $fq2 = "${base}_2.$extension";
        $sample_info{$sample_id} = (-e $fq2) ? { fq1 => $fq1, fq2 => $fq2, paired => 1 } : { fq1 => $fq1, paired => 0 };
    }
}

foreach my $sample_id (keys %sample_info) {
    my $fq1 = $sample_info{$sample_id}{fq1};
    my $fq2 = $sample_info{$sample_id}{fq2} // "";
    my $paired = $sample_info{$sample_id}{paired};

    my $pid = fork();
    if (!defined $pid) {
        die "Cannot fork: $!\n";
    } elsif ($pid == 0) {
        process_sample($sample_id, $fq1, $fq2, $paired);
        exit(0);
    }
}
while (wait() != -1) {}

sub process_sample {
    my ($sample_id, $fq1, $fq2, $paired) = @_;
    my $sample_dir = "$outputdir/$sample_id";
    make_path($outputdir) unless -d $outputdir;
    make_path($sample_dir) unless -d $sample_dir;

    my $trim_cmd = $paired
        ? "env trim_galore --paired --fastqc --cores $threads -o $sample_dir $fq1 $fq2"
        : "env trim_galore --fastqc --cores $threads -o $sample_dir $fq1";
    system($trim_cmd);

    my $accepted_sam = "$sample_dir/accepted_hits.sam";
    my $accepted_sorted_bam = "$sample_dir/accepted_hits.sorted.bam";
    my $accepted_unique_bam = "$sample_dir/accepted_hits.sorted.unique.bam";
    my $metrics_file = "$sample_dir/${sample_id}.metricsFile";
    my $bigwig_output = "$sample_dir/${sample_id}.bw";

    system("echo '======== $sample_id Bowtie2 start ========' >> $bowtie_log_all");
    my $bowtie_cmd = $paired
        ? "bowtie2 -x $indexdir -p $threads -t -q -N 1 -L 25 --no-mixed --no-discordant --rg-id $sample_id --rg SM:$sample_id -1 $fq1 -2 $fq2 -S $accepted_sam >> $bowtie_log_all 2>&1"
        : "bowtie2 -x $indexdir -p $threads -t -q -N 1 -L 25 --no-mixed --no-discordant --rg-id $sample_id --rg SM:$sample_id -U $fq1 -S $accepted_sam >> $bowtie_log_all 2>&1";
    system($bowtie_cmd);
    system("echo '======== $sample_id Bowtie2 end ========' >> $bowtie_log_all");

    system("samtools view -bS $accepted_sam -o $accepted_sorted_bam && samtools sort -@ $threads -o $accepted_sorted_bam $accepted_sorted_bam");
    system("java -Xmx15g -jar $picarddir/picard.jar MarkDuplicates I=$accepted_sorted_bam O=$accepted_unique_bam METRICS_FILE=$metrics_file REMOVE_DUPLICATES=true");
    system("samtools index $accepted_unique_bam");
    system("bamCoverage -b $accepted_unique_bam -of bigwig --binSize 5 --ignoreDuplicates --normalizeUsing BPM --numberOfProcessors $threads -o $bigwig_output");

    my @files_to_delete = (
        "$sample_dir/accepted_hits.sam",
        "$sample_dir/accepted_hits.sorted.bam",
        "$sample_dir/accepted_hits.sorted.bam.bai"
    );
    foreach my $file (@files_to_delete) {
        unlink $file if -e $file;
    }
}


# nohup perl ChIPseq_bam_bwfile.pl --inputdir /home/yxiaobo/RCProj/ChIPseq/database/TF-binding/GSE85632_text/fastqfile --outputdir /home/yxiaobo/RCProj/ChIPseq/database/TF-binding/GSE85632_text --indexdir /home/public_software_annotation/genomeanno/mousebowtie2Index/GRCm38 --picarddir /home/public_software_annotation/software/picard-2.18.2 --threads 8 &
