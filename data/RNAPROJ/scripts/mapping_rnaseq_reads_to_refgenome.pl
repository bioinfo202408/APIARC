#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);

my ($inputdir, $outputdir, $indexdir, $picarddir, $trim_galore_dir, $threads, $mRNA_index);
GetOptions(
    "inputdir=s"        => \$inputdir,
    "outputdir=s"       => \$outputdir,
    "indexdir=s"        => \$indexdir,
    "picarddir=s"       => \$picarddir,
    "trim_galore_dir=s" => \$trim_galore_dir,
    "threads=i"         => \$threads,
    "mRNA_index=s"      => \$mRNA_index,
);

my @supported_formats = ("*fastq", "*fq", "*fastq.gz", "*fq.gz");
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

sub find_trimmed_files {
    my ($dir, $sample_id, $paired) = @_;
    my ($fq1_trimmed, $fq2_trimmed);
    if ($paired) {
        ($fq1_trimmed) = glob("$dir/${sample_id}_1_val_1.*");
        ($fq2_trimmed) = glob("$dir/${sample_id}_2_val_2.*");
        unless ($fq1_trimmed && $fq2_trimmed) {
            die "Error: Unable to find both _val_1 and _val_2 files for sample $sample_id in directory $dir\n";
        }
    } else {
        ($fq1_trimmed) = glob("$dir/${sample_id}_trimmed.*");
        unless ($fq1_trimmed) {
            die "Error: Unable to find _trimmed file for sample $sample_id in directory $dir\n";
        }
    }
    return ($fq1_trimmed, $fq2_trimmed);
}

foreach my $sample_id (keys %sample_info) {
    my $fq1 = $sample_info{$sample_id}{fq1};
    my $fq2 = $sample_info{$sample_id}{fq2} // "";
    my $paired = $sample_info{$sample_id}{paired};
    my $sample_dir = "$outputdir/$sample_id";
    make_path($sample_dir);

    my $pid = fork();
    if ($pid == 0) {
        my $trim_cmd = $paired
            ? "trim_galore --paired --fastqc --cores $threads -o $sample_dir/$trim_galore_dir $fq1 $fq2"
            : "trim_galore --fastqc --cores $threads -o $sample_dir/$trim_galore_dir $fq1";
        system($trim_cmd);

        my ($fq1_trimmed, $fq2_trimmed) = find_trimmed_files("$sample_dir/$trim_galore_dir", $sample_id, $paired);

        my $sample_hisat2_log = "$sample_dir/${sample_id}.hisat2.log";
        my $hisat2_result_log = "$outputdir/hisat2_result.log";
        my $hisat2_cmd = $paired
            ? "hisat2 -x $indexdir -p $threads --dta --rg-id $sample_id --rg SM:$sample_id -1 $fq1_trimmed -2 $fq2_trimmed -S $sample_dir/$sample_id.sam 2> $sample_hisat2_log"
            : "hisat2 -x $indexdir -p $threads --dta --rg-id $sample_id --rg SM:$sample_id -U $fq1_trimmed -S $sample_dir/$sample_id.sam 2> $sample_hisat2_log";
        system($hisat2_cmd);

        if (-e $sample_hisat2_log) {
            open(my $IN, '<', $sample_hisat2_log);
            open(my $OUT, '>>', $hisat2_result_log);
            print $OUT "========== $sample_id ==========\n", <$IN>, "\n";
            close $IN;
            close $OUT;
            unlink $sample_hisat2_log;
        }

        my $filtered_sam = "$sample_dir/$sample_id.filtered.sam";
        system("grep -v -E -w 'NH:i:[2-9]|NH:i:1[0-9]|NH:i:20' $sample_dir/$sample_id.sam > $filtered_sam");

        my $bam_file = "$sample_dir/$sample_id.bam";
        system("samtools view -bS $filtered_sam | samtools sort -@ $threads -o $bam_file");
        system("samtools index $bam_file");

        my $dedup_cmd = qq{java -Xmx15g -jar "$picarddir/picard.jar" MarkDuplicates I=$bam_file O=$sample_dir/$sample_id.dedup.bam METRICS_FILE=$sample_dir/$sample_id.metricsFile REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate};
        system($dedup_cmd);
        system("samtools index $sample_dir/$sample_id.dedup.bam");

        my $mRNA_dir = "$outputdir/mRNA";
        make_path("$mRNA_dir/$sample_id");
        my $stringtie_cmd = "stringtie -p $threads -e -B -G $mRNA_index -A $mRNA_dir/$sample_id/gene_abundance.txt -o $mRNA_dir/$sample_id/transcripts.gtf $sample_dir/$sample_id.dedup.bam";
        system($stringtie_cmd);

        exit(0);
    }
}

while (wait() != -1) {}

# nohup perl mapping_rnaseq_reads_to_refgenome.pl --inputdir /home/yxiaobo/RCProj/RNAseq/database/TF-binding/GSE85632_pipline/fastqfile --outputdir /home/yxiaobo/RCProj/RNAseq/database/TF-binding/GSE85632_pipline --indexdir /home/public_software_annotation/genomeanno/mousehisat2Index/GRCm38 --picarddir /home/public_software_annotation/software/picard-2.18.2  --trim_galore_dir trim_galorefile --mRNA_index /home/public_software_annotation/genomeanno/gencode.vM25.mRNA.annotation.gtf --threads 5 &