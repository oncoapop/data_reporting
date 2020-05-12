# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5 
use strict;
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

my $confFile = "PATH TO YOUR REGISTRY FILE UNDER biomart-perl/conf/. For Biomart Central Registry navigate to
						http://www.biomart.org/biomart/martservice?type=registry";
#
# NB: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
#

my $action='cached';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

		
	$query->setDataset("hsapiens_gene_ensembl");
	$query->addFilter("refseq_mrna", 
["NM_001025088","NM_021938","NM_183425","NM_033010","NM_001130141","NM_032281","NM_014912","NM_007185","NM_005850","NM_031492","NM_005156","NM_021190","NM_006275","NM_005626","NM_018028","NM_006328","NM_018282","NM_001004317","NM_005520","NM_004966","NM_017697","NM_021938","NM_005839","NM_016333","NM_006445","NM_176800"]);
	$query->addAttribute("ensembl_transcript_id");
	$query->addAttribute("ensembl_peptide_id");

$query->formatter("TSV");

my $query_runner = BioMart::QueryRunner->new();
############################## GET COUNT ############################
# $query->count(1);
# $query_runner->execute($query);
# print $query_runner->getCount();
#####################################################################


############################## GET RESULTS ##########################
# to obtain unique rows only
# $query_runner->uniqueRowsOnly(1);

$query_runner->execute($query);
$query_runner->printHeader();
$query_runner->printResults();
$query_runner->printFooter();
#####################################################################
