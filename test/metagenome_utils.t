use Test::Most;
use Data::Dumper::Concise;

use_ok 'MetagenomeUtils';

subtest 'testing build_prodigal_params' => sub {

    my %defaults = (
        "-c"  => '',
        "-f"  => 'sco',
        "-g"  => 11,
        "-p"  => 'meta',
        "-q"  => '',
    );

    throws_ok {
        MetagenomeUtils::build_prodigal_params();
    }   qr!An input FASTA/Genbank file is required for Prodigal to run!,
        'build_prodigal_params dies without args';

    ##

    my @test_args = ({
        input   => [ 'blob.fasta' ],
        output  => { '-i' => 'blob.fasta' },
    },{

        input   => [ 'blob.fasta', 'trans.file' ],
        output  => { '-i' => 'blob.fasta', '-a' => 'trans.file' },
    });

    for my $test ( @test_args ) {

        cmp_deeply
            MetagenomeUtils::build_prodigal_params( @{ $test->{ input } } ),
            { %defaults, %{ $test->{ output } } },
            'appropriate command line params with args ' . Dumper $test->{ input }
            or diag explain {
                got     =>
                expect  => { %defaults, %{ $test->{ output } } },
            };
    }


};

subtest 'run_prodigal' => sub {

    throws_ok {
        MetagenomeUtils::run_prodigal()
    } qr/No parameters supplied to run_prodigal/,
        'run_prodigal dies without params';

    throws_ok {
        MetagenomeUtils::run_prodigal( {} )
    } qr/No parameters supplied to run_prodigal/,
        'run_prodigal dies with an empty hashref';

    # override system locally to check what params we have been given
    {
        no warnings 'redefine';
        local *MetagenomeUtils::_syscall = sub {
            cmp_deeply
                \@_,
                [ '/kb/runtime/prodigal', '-i', 'fasta.file' ],
                'Got the expected input array';
            # mimic a successful response
            0;
        };

        lives_ok {
            MetagenomeUtils::run_prodigal({ '-i' => 'fasta.file' });
        } 'run prodigal runs OK with appropriate args';

    }

    {
        no warnings 'redefine';
        local *MetagenomeUtils::_syscall = sub {
            # mimic death
            'Oh no! I died';
        };

        throws_ok {
            MetagenomeUtils::run_prodigal({ '-i' => 'fasta.file' });
        } qr/Prodigal run failed: Oh no! I died/, 'prodigal dies appropriately';

    }

};


my $file_not_found = '/this_file/does/not/exist';

subtest 'parse_protein_translation_file' => sub {

    throws_ok {
        MetagenomeUtils::parse_protein_translation_file( $file_not_found );
    } qr/Could not read-open $file_not_found:/,
        'file not found error';

    # create a file with an unparsable record
#     throws_ok {
#         MetagenomeUtils::parse_protein_translation_file( '/path/to/dodgy/file' );
#     } qr/Could not parse record:/,
#         'could not parse file with incorrect record';

    my $expected = {
        # whatever the expected data structure is
    };

    # successful parse
#     cmp_deeply
#         MetagenomeUtils::parse_protein_translation_file( '/path/to/good/file' ),
#         $expected,
#         'got the expected data structure from parsing the translation file';

};

subtest 'parse_prodigal_output_file' => sub {

    throws_ok {
        MetagenomeUtils::parse_prodigal_output_file( $file_not_found );
    } qr/Could not read-open $file_not_found:/,
        'file not found error';

#     my $expected = {
#         # whatever the expected data structure is
#     };
#
#     # successful parse
#     cmp_deeply
#         MetagenomeUtils::parse_prodigal_output_file( '/path/to/good/file' ),
#         $expected,
#         'got the expected data structure from parsing the prodigal output file';

};

done_testing;
