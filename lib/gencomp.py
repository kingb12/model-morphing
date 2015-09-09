from biokbase.GenomeComparison.Client import GenomeComparison
gencomp = GenomeComparison(url='https://next.kbase.us/services/genome_comparison/')
gencomp.blast_proteomes({'genome1ws' : '68', 'genome2ws' : '68', 'genome1id' : '2',
                         'genome1id' : '3', 'output_ws' : '68', 'output_id' : 'acetivo-barkeri-comp'})



