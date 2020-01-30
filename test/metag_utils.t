use Test::Most;
use Data::Dumper::Concise;

use_ok 'metag_utils';

subtest 'testing build_prodigal_params' => sub {

    my @default = ( "-c", "-f sco", "-p meta", "-q" );


    throws_ok {
        metag_utils::build_prodigal_params();
    }   qr!An input FASTA/Genbank file is required for Prodigal to run!,
        'build_prodigal_params dies without args';

    ##

    my @test_args = ({
        input   => [ 'this', 'blob.fasta' ],
        output  => [ '-i blob.fasta' ],
    },{

        input   => [ 'this', 'blob.fasta', 'trans.file' ],
        output  => [ '-i blob.fasta', '-a trans.file' ],
    });

    for my $test ( @test_args ) {

        cmp_deeply
            metag_utils::build_prodigal_params( @{ $test->{ input } } ),
            [ @default, @{ $test->{ output } } ],
            'appropriate command line params with args ' . Dumper $test->{ input }
            or diag explain {
                got     =>
                expect  => [ @default, @{ $test->{ output } } ],
            };
    }


};

subtest 'run_prodigal' => sub {

    throws_ok {
        metag_utils::run_prodigal()
    } qr/No parameters supplied to run_prodigal/,
        'run_prodigal dies without params';

};

subtest 'parse_prodigal_results' => sub {

    ok 1;
};

done_testing;