#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Statistics::Descriptive;
use Scalar::Util qw(looks_like_number);

# ================================================== #
#      Name : raw_pheno_outlier_filter                      
#    Author : Liu Peng                               #
#    E-mail : sxliulian2012@hotmail.com              #
#   Version : 0.1 								     #
#      Date : 2021-11-27 14:28:20                    #
# ================================================== #
# dependent package:
# Statistics::Descriptive
# Scalar::Util

die "Usage: perl $0 input.txt [prefix] \n" if @ARGV < 1;
my $prefix = $ARGV[1] ||= "out";

open IN, $ARGV[0] || die $!;
my $header = <IN>;
$header =~ s/[\r\n]//g;

my @traits = split /\t/, $header;

my %data_hash;
my %id_hash;

while(<IN>){
    $_ =~ s/[\r\n]//g;
    next if /^$|^#/;

    my @a = split /\t/, $_;
    for (my $i = 1; $i<@a; $i++){
        push @{$data_hash{$traits[$i]}{$a[0]}}, $a[$i];
    }

    $id_hash{$a[0]}++;
}
close IN;

my %mean;
my %out_good;
my %out_bad;

foreach my $trait (sort keys %data_hash){
    foreach my $id (sort keys %{$data_hash{$trait}}){
        my $ddd = $data_hash{$trait}{$id};
        
        # check data
        my (@pass, @drop);
        foreach my $tmp1 (@$ddd){
            if (! looks_like_number($tmp1)){
                push @drop, $tmp1;
            }elsif ($tmp1 == 0) {
                push @drop, $tmp1;
            }else{
                push @pass, $tmp1;
            }
        }
        my @good;
        # check outlier
        if (scalar(@pass) > 3){
            # my ($g, $b, $m) = &filter_by_3u(@pass);
            my ($g, $b, $m) = &filter_by_zscore(\@pass, 1);
            

            $mean{$id}{$trait} = $m;
            push @good, @$g;
            push @drop, @$b;

        }elsif(scalar(@pass) > 0){
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@pass);
            my $mean = $stat->mean();
            
            $mean{$id}{$trait} = $mean;
            push @good, @pass;

        }else{
            next;
        }
        $out_good{$id}{$trait} = join("|", @good);
        $out_bad{$id}{$trait} = join("|", @drop);

    }
}

# print Dumper %mean;
## output mean
open O1, "> $prefix.mean.txt" || die $!;
print O1 join("\t", "Taxa", sort @traits[1..$#traits])."\n";
# sort by id
my @arr = keys %mean;
my @p_sort = map{$_->[0]}
    sort {$a->[1] <=> $b->[1] }    
        map{ my ($idx) = $_ =~ /[IBMSCL](\d+)/; [$_, $idx]; }@arr;


# foreach my $ibm (sort {$a<=>$b} keys %mean){
foreach my $ibm (@p_sort){
    my @out;
    foreach my $ttt (sort @traits[1..$#traits]){
        my $mmm = exists $mean{$ibm}{$ttt} ? sprintf("%.2f", $mean{$ibm}{$ttt}) : "NA";
        push @out, $mmm;
    }
    print O1 join("\t", $ibm, @out)."\n";
}
close O1;

## output good data
open O2, "> $prefix.good.txt" || die $!;
print O2 join("\t", "Taxa", sort @traits[1..$#traits])."\n";

foreach my $ibm (@p_sort){
    my @out;
    foreach my $ttt (sort @traits[1..$#traits]){
        my $mmm = exists $out_good{$ibm}{$ttt} ? $out_good{$ibm}{$ttt} : "NA";
        push @out, $mmm;
    }
    print O2 join("\t", $ibm, @out)."\n";
}

close O2;

# output bad data
open O3, "> $prefix.bad.txt" || die $!;
print O3 join("\t", "Taxa", sort @traits[1..$#traits])."\n";

foreach my $ibm (@p_sort){
    my @out;
    foreach my $ttt (sort @traits[1..$#traits]){
        my $mmm = exists $out_bad{$ibm}{$ttt} ? $out_bad{$ibm}{$ttt} : "NA";
        push @out, $mmm;
    }
    print O3 join("\t", $ibm, @out)."\n";
}

close O3;

###
sub filter_by_3u {
    my @ddd = @_;

    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@ddd);
    my $Q1 = $stat->quantile(1);
    my $Q3 = $stat->quantile(3);
    my $e_value_min = $Q1 - 1.5 * ($Q3 - $Q1);
    my $e_value_max = $Q3 + 1.5 * ($Q3 - $Q1);

    my @good;
    my @bad;

    foreach my $value (@ddd){
        if ($value >= $e_value_min && $value <= $e_value_max ){
            push @good, $value;
        }else{
            push @bad, $value;
        }
    }
    ## calc mean
    my $tmp = Statistics::Descriptive::Full->new();
    $tmp->add_data(@good);
    my $m = $tmp->mean();

    return (\@good, \@bad), $m;
}

sub filter_by_zscore {
    # filter by zscore
    my ($data_arr, $cutoff) = @_;
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@$data_arr);
    my $mean = $stat->mean();
    my $standard_deviation = $stat->standard_deviation();
    
    my @good;
    my @bad;

    foreach my $value (@$data_arr){
        my $zscore = ($value - $mean) / $standard_deviation;
        if ( abs($zscore) <= $cutoff ) {
            push @good, $value;

        }else{
            push @bad, $value;
        }
    }
    ## calc mean
    my $tmp = Statistics::Descriptive::Full->new();
    $tmp->add_data(@good);
    my $m = $tmp->mean();

    return (\@good, \@bad, $m);
}

__END__

data format:
Taxa trait1 trait2 ...
O1_rep1 1   2.3
O1_rep2 2.3 2.3
O1_rep3 .. .. 
O1_rep4 .. ..
O2_rep1 .. ..
O2_rep2 .. ..
O2_rep3 .. ..