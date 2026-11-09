#!/usr/bin/env perl
=head1 NAME

generate_tree.pl - Generate synthetic phylogenetic trees strictly using BioPerl Bio::Tree::RandomFactory

=head1 SYNOPSIS

  perl src/generate_tree.pl --taxa 32 --scale 1.0 --seed 42 --outtree tree.nwk

=head1 DESCRIPTION

Generates synthetic phylogenetic trees according to Matsui & Iwasaki (2020, Systematic Biology):
- Topology & Branch Lengths: BioPerl Bio::Tree::RandomFactory using the backward Yule process (rand_yule_c_tree).
- Model: Constant birth rate lambda = 1.0, branch lengths sampled from 1 - ln(u * (e - 1) + 1).
- Scaling: Multiplied by global evolutionary distance scale D (--scale).
- Output: Standard Newick format without internal node labels.

=cut

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use Bio::Tree::RandomFactory;
use Bio::TreeIO;

my $taxa      = 32;
my $scale     = 1.0;
my $seed      = undef;
my $model     = "paper_yule";
my $rate_sd   = 0.0;
my $lba_ratio = 1.0;
my $outtree   = undef;
my $help      = 0;

GetOptions(
    "taxa=i"      => \$taxa,
    "scale=f"     => \$scale,
    "seed=i"      => \$seed,
    "model=s"     => \$model,
    "rate_sd=f"   => \$rate_sd,
    "lba_ratio=f" => \$lba_ratio,
    "outtree=s"   => \$outtree,
    "help|h"      => \$help,
) or die("Error in command line arguments\n");

if ($help || !$outtree) {
    print "Usage: $0 --taxa <N> --scale <D> --seed <S> --outtree <file.nwk>\n";
    exit($help ? 0 : 1);
}

die "Taxa count must be at least 2, got $taxa\n" if $taxa < 2;
die "Scale must be strictly positive, got $scale\n" if $scale <= 0.0;

# Seed random number generator if provided
if (defined $seed) {
    srand($seed);
}

# 1. Initialize taxon names: T1, T2, ..., TN
my @taxon_names = map { "T$_" } (1..$taxa);

# 2. Generate tree using BioPerl Bio::Tree::RandomFactory (backward Yule process)
my $factory = Bio::Tree::RandomFactory->new(
    -taxa     => \@taxon_names,
    -randtype => "yule",
);

my $tree = $factory->next_tree;

# Helper for Gaussian random number (Box-Muller transform)
sub gaussian_rand {
    my ($mean, $sd) = @_;
    return $mean if $sd <= 0.0;
    my $u1 = rand();
    while ($u1 <= 1e-15) { $u1 = rand(); }
    my $u2 = rand();
    my $z = sqrt(-2.0 * log($u1)) * cos(2.0 * 3.14159265358979323846 * $u2);
    return $mean + $z * $sd;
}

# 3. Apply scaling, rate heterogeneity, LBA, and clear internal node IDs
my $min_length = 1e-6;

for my $node ($tree->get_nodes) {
    if (!$node->is_Leaf) {
        $node->id("");
    }

    my $bl = $node->branch_length;
    if (defined $bl) {
        my $scaled_bl = $bl * $scale;

        if ($rate_sd > 0.0) {
            my $perturbation = gaussian_rand(0.0, $rate_sd * $scale);
            $scaled_bl = $scaled_bl + $perturbation;
        }

        if ($lba_ratio > 1.0 && $node->is_Leaf) {
            my $id = $node->id || "";
            if ($id eq "T1" || $id eq "T$taxa") {
                $scaled_bl *= $lba_ratio;
            }
        }

        if ($scaled_bl < $min_length) {
            $scaled_bl = $min_length;
        }
        $node->branch_length($scaled_bl);
    }
}

# 4. Serialize to Newick format using Bio::TreeIO
my $newick_str = "";
open(my $io, ">", \$newick_str) or die "Cannot open in-memory scalar buffer: $!\n";
my $treeio = Bio::TreeIO->new(-format => "newick", -fh => $io);
$treeio->write_tree($tree);
close($io);

$newick_str =~ s/^\s+|\s+$//g;
if ($newick_str !~ /;$/) {
    $newick_str .= ";";
}

my $out_dir = dirname($outtree);
if ($out_dir && ! -d $out_dir) {
    make_path($out_dir);
}

open(my $fh, ">", $outtree) or die "Cannot write to $outtree: $!\n";
print $fh "$newick_str\n";
close($fh);

