import sys
import os
import tempfile
import optparse
import subprocess

from json import loads, dumps

DEFAULT_DATA_TABLE_NAME = "ensembl_cache"

def fetch_ensembl_cache(data_manager_dict, ensembl_vep_path, species, assembly, ensembl_version, params, target_directory, data_table_name):

    perl_executable = '/'.join([ensembl_vep_path, "INSTALL.pl"])

    commands = [ "perl" , perl_executable, "--AUTO", "cf", "--SPECIES", species, "--ASSEMBLY", assembly, "--CACHEDIR", target_directory ]

    print commands

    proc = subprocess.Popen( args=commands, shell=False, cwd=target_directory)
    return_code = proc.wait()
    if return_code:
        print >> sys.stderr, "Error fetching ensembl cache"
        sys.exit( return_code )
    
    data_table_entry = dict(value='.'.join([species,assembly,ensembl_version]), name='.'.join([species,assembly,ensembl_version]), species=species, assembly=assembly, version=ensembl_version, path=target_directory )
    
    _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry )

def _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][ data_table_name ] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][ data_table_name ].append( data_table_entry )
    return data_manager_dict
    
def main():
    
    print >> sys.stderr, "Starting" 
    
    parser = optparse.OptionParser()
    parser.add_option( '-i', '--organism_string', dest='input', action='store', type="string")
    parser.add_option( '-e', '--ensembl_vep_path', dest='ensembl_vep_path', action='store', type='string')
    (options, args) = parser.parse_args()
    
    filename = args[0]

    print >> sys.stderr, "What"
    print >> sys.stderr, "When"

    ensembl_vep_path = options.ensembl_vep_path
    species = options.input.split("-")[0]
    assembly = options.input.split("-")[1]
    ensembl_version = 84

    print >> sys.stderr, "hello"

    params = loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    data_manager_dict = {}
    
    #Fetch the Ensembl Cache
    fetch_ensembl_cache( data_manager_dict, ensembl_vep_path, species, assembly, ensembl_version, params, target_directory, data_table_name=DEFAULT_DATA_TABLE_NAME )
       
 
    #save info to json file
    open( filename, 'wb' ).write( dumps( data_manager_dict ) )

if __name__ == "__main__": main()        
