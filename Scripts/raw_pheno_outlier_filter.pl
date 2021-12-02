#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Statistics::Descriptive;
use Scalar::Util qw(looks_like_number);
use Getopt::Long;
use File::Basename;

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

my ($input, $output_prefix, $method, $cutoff, $good, $bad, $help);
GetOptions(
    'i|input=s'     =>  \$input,
    'o|output:s'    =>  \$output_prefix,
    'm|method:s'    =>  \$method,
    'c|cutoff:f'    =>  \$cutoff,
    'g|good!'       =>  \$good,
    'b|bad!'        =>  \$bad,
    'h|help'        =>  \$help,
) or die $!;

die &usage() if (defined $help or !defined $input);
$output_prefix ||= "./output";  # prefix of output with path.
$method ||= 'zscore';           # zscore or 3u
$cutoff ||= '1.5';              # cutoff for filter.

my ($prefix, $outdir, $suffix) = fileparse($output_prefix);
if (! -d $outdir){
    mkdir $outdir;
}
my $output_mean_file = $outdir . $prefix . "_mean.txt";
print STDOUT "Output mean file: $output_mean_file \n";
## begin
my ($data_hash, $row_ids, $colum_ids) = &data_to_hash($input);

## filter
my (%out_mean, %out_good, %out_bad);
foreach my $col (sort keys %$data_hash){
    foreach my $row (sort keys %{$data_hash->{$col}}){
        my @values = @{$data_hash->{$col}->{$row}};

        # check values
        my (@pass_val, @good_val, @drop_val);

        foreach my $val (@values){
            if (looks_like_number($val) && $val != 0){
                push @pass_val, $val;
            }else{
                push @drop_val, $val;
            }
        }

        if (scalar(@pass_val) >= 3){
            # filter
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@pass_val);
            my $mean = $stat->mean();

            if ($method eq '3u'){
                my $Q1 = $stat->quantile(1);
                my $Q3 = $stat->quantile(3);
                my $e_value_min = $Q1 - $cutoff * ($Q3 - $Q1);
                my $e_value_max = $Q3 + $cutoff * ($Q3 - $Q1);

                foreach my $val (@pass_val){
                    if ($val >= $e_value_min && $val <= $e_value_max){
                        push @good_val, $val;
                    }else{
                        push @drop_val, $val;
                    }
                }
            }elsif($method eq 'zscore'){
                my $standard_deviation = $stat->standard_deviation() ? $stat->standard_deviation() : 0;

                if ($standard_deviation > 0){
                    foreach my $val (@pass_val){
                        my $zscore = ($val - $mean) / $standard_deviation;
                        if (abs($zscore) <= $cutoff ) {
                            push @good_val, $val;
                        }else{
                            push @drop_val, $val;
                        }
                    }
                }else{
                    push @good_val, @pass_val;
                }
            }else{
                print STDERR "ERROR!\n";
                exit 1;
            }

            $out_mean{$row}{$col} = &average(@good_val);
            $out_good{$row}{$col} = join("|", @good_val);

        }elsif (scalar(@pass_val) > 0) {
            push @good_val, @pass_val;
            $out_mean{$row}{$col} = &average(@good_val);
            $out_good{$row}{$col} = join("|", @good_val);
        }else{
            $out_mean{$row}{$col} = "NA";
            $out_good{$row}{$col} = "-";
            $out_bad{$row}{$col}  = join("|", @drop_val);
        }
        # $out_mean{$row}{$col} = &average(@good_val);
        # $out_good{$row}{$col} = join("|", @good_val);
        $out_bad{$row}{$col}  = scalar(@drop_val) > 0 ? join("|", @drop_val) : "-"; 
    }
}


# print Dumper %out_mean;
# sort by row id
my @p_sort = map{$_->[0]}
    sort {$a->[1] <=> $b->[1] }    
        map{ my ($idx) = $_ =~ /[IBMSCL](\d+)/; [$_, $idx]; } keys %$row_ids;

# output mean values file
open O1, "> $output_mean_file" || die "Can't write to file: $output_mean_file \n";
print O1 join("\t", "Taxa", sort keys %$colum_ids)."\n";

foreach my $row_id (@p_sort){
    my @out;
    foreach my $col_id (sort keys %$colum_ids){

        # my $mmm = exists $out_mean{$row_id}{$col_id} ? sprintf("%.2f", $out_mean{$row_id}{$col_id}) : "NA";
        my $mmm = exists $out_mean{$row_id}{$col_id} ? $out_mean{$row_id}{$col_id} : "NA";

        push @out, $mmm;
    }
    print O1 join("\t", $row_id, @out)."\n";
}
close O1;

if (defined $good){
    my $output_good_file = $outdir . $prefix . "_good.txt";
    print STDOUT "Output good values file: $output_good_file \n";

    # output good data
    open O2, "> $output_good_file" || die "Can't write to the file: $output_good_file\n";
    print O2 join("\t", "Taxa", sort keys %$colum_ids)."\n";
    foreach my $row_id (@p_sort){
        my @out;
        foreach my $col_id (sort keys %$colum_ids){
            my $mmm = exists $out_good{$row_id}{$col_id} ? $out_good{$row_id}{$col_id} : "NA";
            push @out, $mmm;
        }
        print O2 join("\t", $row_id, @out)."\n";
    }
    close O2;
}

if (defined $bad){
    my $output_bad_file = $outdir . $prefix . "_bad.txt";
    print STDOUT "Output good values file: $output_bad_file \n";

    # output good data
    open O3, "> $output_bad_file" || die "Can't write to the file: $output_bad_file\n";
    print O3 join("\t", "Taxa", sort keys %$colum_ids)."\n";
    foreach my $row_id (@p_sort){
        my @out;
        foreach my $col_id (sort keys %$colum_ids){
            my $mmm = exists $out_bad{$row_id}{$col_id} ? $out_bad{$row_id}{$col_id} : "NA";
            push @out, $mmm;
        }
        print O3 join("\t", $row_id, @out)."\n";
    }
    close O3;
}

## ......................................................................... ##
sub data_to_hash {
    # read raw data to hash
    my $infile = shift @_;;

    open IN, "< $infile" || die "Can't open file: $infile \n";
    my $header = <IN>; $header =~ s/[\r\n]//g;
    my @colums = split /\t/, $header;
    my (%data_hash, %row_id_hash, %colum_id_hash);

    # record uniq col names
    foreach (@colums[1..$#colums]){$colum_id_hash{$_} = 1};

    while(<IN>){
        $_ =~ s/[\r\n]//g; next if /^$|^#/;
        my @a = split /\t/, $_;

        for (my $i = 1; $i<@a; $i++){
            push @{$data_hash{$colums[$i]}{$a[0]}}, $a[$i];
        }

        $row_id_hash{$a[0]}++;
    }
    close IN;

    return (\%data_hash, \%row_id_hash, \%colum_id_hash);
}

sub check {
    # check data
    my @data = @_;
    my %out;

    foreach my $v (@data){
        if (!looks_like_number($v)){
            push @{$out{'drop'}}, $v;
            next;
        }elsif ($v == 0) {
            push @{$out{'drop'}}, $v;
            next;
        }else{
            push @{$out{'pass'}}, $v;
            next;
        }
    }
    return(\%out);
}

sub data_filter {
    # set filter method
    my ($arr_values, $filter_method, $filter_cutoff) = @_;
    # print Dumper $arr_values;
    my $checked_data = &check($arr_values);
    my @pass_val = $checked_data->{'pass'};
    my @drop_val = $checked_data->{'drop'};
    
    my (%out, @good_val, @tmp);
    if (scalar (@pass_val) >= 3){
        if ($method eq "zscore"){
            (@good_val, @tmp) = &filter_by_zscore(\@pass_val, $filter_cutoff);

        }elsif($method eq "3u"){
            (@good_val, @tmp) = &filter_by_3u(\@pass_val, $filter_cutoff);

        }else{
            print STDERR "Please check method params, only 'zscore' or '3u' can be identify\n";
            exit 1; 
        }
        push @drop_val, @tmp;

    }elsif(scalar (@pass_val) > 0){
        push @good_val, @pass_val;

    }else{
        push @good_val, "NA";
    }

    if (scalar (@drop_val) == 0 || undef @drop_val){
        push @drop_val, "NA";
    }

    push @{$out{'good'}}, @good_val;
    push @{$out{'drop'}}, @drop_val;

    return(\%out);
}

sub average {
    # calculate the mean value of data
    use Statistics::Descriptive;

    my @data = @_;
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@data);
    my $mean = $stat->mean();

    $mean = sprintf("%.2f", $mean);

    return($mean);
}

# sub filter_by_3u {
#     # filter by 3u
#     use Statistics::Descriptive;

#     my ($data, $cutoff) = @_;
#     my $stat = Statistics::Descriptive::Full->new();
#     $stat->add_data(@$data);
#     my $Q1 = $stat->quantile(1);
#     my $Q3 = $stat->quantile(3);
#     my $e_value_min = $Q1 - $cutoff * ($Q3 - $Q1);
#     my $e_value_max = $Q3 + $cutoff * ($Q3 - $Q1);

#     my %out;

#     my (@good_val, @drop_val);

#     foreach my $value (@$data){
#         if ($value >= $e_value_min && $value <= $e_value_max ){
#             push @good_val, $value;
#         }else{
#             push @drop_val, $value;
#         }
#     }
#     push @{$out{'good'}}, @good_val;
#     push @{$out{'drop'}}, @drop_val;

#     return (\%out);
# }

# sub filter_by_zscore {
#     # filter by zscore
#     use Statistics::Descriptive;

#     my ($data, $cutoff) = @_;
#     my $stat = Statistics::Descriptive::Full->new();
#     $stat->add_data(@$data);
#     my $mean = $stat->mean();
#     my $standard_deviation = $stat->standard_deviation();
#     my (%out, @good_val, @drop_val);

#     # print Dumper $standard_deviation;

#     # if ($standard_deviation > 0){
#     #     foreach my $value (@$data){
#     #         my $zscore = ($value - $mean) / $standard_deviation;
#     #         if (abs($zscore) <= $cutoff ) {
#     #             push @good_val, $value;
#     #         }else{
#     #             push @drop_val, $value;
#     #         }
#     #     }
#     # }else{
#     #     push @good_val, @$data;
        
#     # }

#     # push @{$out{'good'}}, @good_val;
#     # push @{$out{'drop'}}, @drop_val;

#     # return(\%out);
# }

sub usage {
    print STDERR <<EOF;
  Name: 
    $0
  Description: 
    Filter outlier data from raw data.
  Usage: 
    perl $0 [options] -i <input.txt> 
  Options:
    -h|--help   print this help information.
    -i|--input  input raw data.
    -m|--method method for data filter, 3u or zscore. [default: zscore]
    -c|--cutoff cutoff for filter.  [default: 1.5]
    -o|--output prefix of output.   [default: output]
    -g|--good   output of good data.
    -b|--bad    output of bad data.
  Example:
    perl $0 -m 3u -c 1.5 -i input.txt -o /path/to/output -g 
  or
    perl $0 -m zscore -c 1 -i input.txt  -o /path/to/output -g -b 
EOF
    exit 1;
}


__END__

data format:
Taxa trait1 trait2 ...
O1_rep1 1   2.3
O1_rep2 2.3 2.3
O1_rep3 .. .. 
O2_rep1 .. ..
O2_rep2 .. ..

