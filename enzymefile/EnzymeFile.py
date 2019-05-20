import re
import pprint
import sys

class :
    """
    This class deals with data from the KEGG genes/genome file.

    genome file is a DBGET file containing all the genomes (organisms) entries and its related data.

    At this momment KEGG genome file is a small file (3MB). Because of that this class doesn't provide any
    interative way to read the file. The whole file is load in memory by a dicionary data structure. In a machine
    with 4GB RAM and Intel Core i5 it tooks less than 10 seconds to load the data.

    Attributes:
        organism_entries(dict): Where actual the Genome entries are.
        tax_ids_and_its_organisms(list): List of taxonomy ids and its organisms.
        organism_code_and_its_tax_ids(dict): Organism KEGG codes and its taxonomy ids.
        taxonomy_id_and_its_organism_data(list): List of genome data based on its taxonomy id.
        internal_kegg_ids_and_its_organisms(dict): Internal KEGG id and its organism data.
     
    """

    def __init__(self, reader):
        self.dbgetr = reader 

        # Those are important variables because it's supposed to query them a lot of times in different contexts when the actual
        # mean of this class (be part of a KEGG importing process and its tons of text files) is fulfilled. Isolated, those variables seems to be a bad idea.
        self.organism_entries = {}
        self.tax_ids_and_its_organisms = {}
        self.organism_code_and_its_tax_ids = {}
        self.taxonomy_id_and_its_organism_data = {}
        self.internal_kegg_ids_and_its_organisms = {}

    def generate_genome_data( self ):
        """
        Generate a dictionary containing genome file data.

        This method is selective, in other words, despite the 'genome' file has a great amount of information,
        this method only append some specific (and most useful) data.

        But feel free to improve the result dictionary with more fields from the 'genome' file.

        Returns:
            (dict): Genome file data in a dictionary format.
        """

        genome_entries = self.dbgetr.all_entries()

        re_organism_internal_id = re.compile('^(T[0-9]{1,})\s(.*)$')

        for genome_entry in genome_entries:

            # NAME is a typical key in the 'genome' file.
            # And Draft Genomes doesn't have a key NAME.
            # So, we're using 'NAME' as a criteria to tell if it's a valid genome entry.
            if 'NAME' in genome_entry:

                # Get the three letters Kegg organism code.
                organism_code = genome_entry['NAME'][0].split(',')
                organism_code = organism_code[0]
                organism_code = organism_code.replace(' ','')

                # Get the taxonomy NCBI id.
                taxonomy_id = genome_entry['TAXONOMY'][0]
                taxonomy_id = taxonomy_id.split(':')
                taxonomy_id = taxonomy_id[1]

                # Get the definition name gave from KEGG
                definition_name = genome_entry['DEFINITION'][0]
 
                # Get the internal KEGG id for the organism (that's not the three letter code, but something like 'T01445').
                organism_internal_id = genome_entry['ENTRY'][0]
                organism_internal_id_search = re_organism_internal_id.search( organism_internal_id )
                organism_internal_id = organism_internal_id_search.group(1)

                genome_status = organism_internal_id_search.group(2)
                genome_status = re.sub('^\ {1,}', '', genome_status )

                status = re.compile('^Complete')
                status_result = status.search( genome_status )
                if status_result:
                    genome_status = True
                else:
                    genome_status = False

                data = { 'organism_code': organism_code, 'kegg_definition_name': definition_name, 'taxonomy_id': taxonomy_id, 'kegg_organism_id': organism_internal_id, 'complete_genome': genome_status } 

                self.internal_kegg_ids_and_its_organisms[ organism_internal_id ] = organism_code 
                self.organism_entries[ organism_code ]                      = data 

                # TODO: taxonomy id is not uniq!!!!
                if not taxonomy_id in self.taxonomy_id_and_its_organism_data:
                    self.taxonomy_id_and_its_organism_data[ taxonomy_id ] = []

                self.taxonomy_id_and_its_organism_data[ taxonomy_id ].append( data ) 

                # TODO: taxonomy id is not uniq!!!!
                if not taxonomy_id in self.tax_ids_and_its_organisms:
                    self.tax_ids_and_its_organisms[ taxonomy_id ] = []

                self.tax_ids_and_its_organisms[ taxonomy_id ].append( organism_code )

                # TODO: taxonomy id is not uniq!!!!
                self.organism_code_and_its_tax_ids[ organism_code ] = taxonomy_id  

        return self.organism_entries 


    def load_genome_data( self ):
        """
        Load the Genome file in case it wasn't already loaded.
        """

        if len( self.organism_entries ) == 0:
            self.generate_genome_data()


    def all_internal_kegg_organism_ids( self ):
        """
        Return a dictionary containing all KEGG internal organism ids (Ex: T01001, T01453) by KEGG three (or four) letter organism codes.

        Returns:
            (dict): Example:  { 'xfs': 'T03318', 'xft': 'T00113', 'xfu': 'T03094', 'xtr': 'T01011', 'zmc': 'T03385' }
        """

        self.load_genome_data()

        organism_and_ids = {}

        for org,values in self.organism_entries.items():
            organism_and_ids[ org ] = values['kegg_organism_id']

        return organism_and_ids 


    def organism_code_by_internal_kegg_id( self, kegg_id=None ):
        """
        Returns the three (or four) letters KEGG organism code (Ex: hsa, lma, pbn) by the internal KEGG orgnanism id (Ex: T01001, T01023).

        Args:
            kegg_id(str): KEGG internal organism id (Ex: T01001, T02310).

        Returns:
            (str): KEGG orgnanism code (Ex: hsa, ptb, lma, tcr).
        """

        self.load_genome_data()

        if str(kegg_id) in self.internal_kegg_ids_and_its_organisms:
            return self.internal_kegg_ids_and_its_organisms[ str(kegg_id) ]
        else:
            return None


    def all_genome_data( self ):
        """
        Returns all the Genome file in a dictionary format.
        """

        self.load_genome_data()

        return self.organism_entries


    def genome_data_by_genome_code( self, organism_code=None ):
        """
        Returns the Genome data by the three (or four) letter KEGG organism code (hsa, lma, tcr, etc)

        Args:
            organism_code(str): KEGG organism code (hsa, lma, tcr)

        Returns:
            (dict): Dictionary containing the data from the selected Genome organism.
        """

        self.load_genome_data()

        if organism_code in self.organism_entries:
            return self.organism_entries[ organism_code ]
        else:
            return {}


    def genome_data_by_taxonomy_id( self, taxonomy_id=None ):
        """
        Returns the Genome data by the taxonomy id (same as in the NCBI taxonomy database).

        Args:
            taxonomy_id(str): NCBI taxonomy id (Ex: 9606, 3454) 

        Returns:
            (dict): Dictionary containing the data from the selected Genome organism.
        """

        self.load_genome_data()

        if str(taxonomy_id) in self.taxonomy_id_and_its_organism_data:
            return self.taxonomy_id_and_its_organism_data[ str(taxonomy_id) ]
        else:
            return {}


    def all_organism_codes( self ):
        """
        Returns a list containing all the three (or four) letters KEGG organism code (Ex: hsa, tcr, lma, etc).

        Returns:
            (list): KEGG organism codes.
        """

        self.load_genome_data()

        organism_codes = []

        for genome_entry in self.organism_entries:
            organism_codes.append( genome_entry ) 

        return organism_codes 

    # TODO: TAXONOMY ID IS NOT UNIQUE!! FIX IT!!!
    def all_taxonomy_ids( self ):
        """
        Returns a list containing all the NCBI taxonomy ids found in the Genome file (ONLY THOSE FOUND IN THE GENOME FILE. This is not a method that returns NCBI full database taxonomy data).

        Returns:
            (list): NCBI taxonomy ids. 
        """

        if len( self.organism_entries ) == 0:
            self.generate_genome_data()

        taxonomy_ids = []

        for genome_code,data in self.organism_entries.items():
            taxonomy_ids.append( data['taxonomy_id'] ) 

        return taxonomy_ids 

    # TODO: TAXONOMY ID IS NOT UNIQUE!! FIX IT!!!
    def taxonomy_id_by_organism_code( self, organism_code=None ):
        """
        Returns NCBI taxonomy id associated with KEGG organism code.

        Args:
            organism_code(str): Three or four letter KEGG organism code (Ex: hsa, lma, etc)

        Returns:
            (str): NCBI taxonomy id.
        """

        self.load_genome_data()

        if organism_code in self.organism_code_and_its_tax_ids:
            return self.organism_code_and_its_tax_ids[ organism_code ]
        else:
            return None

    # TODO: TAXONOMY ID IS NOT UNIQUE!! FIX IT!!!
    def organism_code_by_taxonomy_id( self, tax_id=None ):
        """
        Returns KEGG organism code associated with NCBI taxonomy id.

        Args:
            tax_id(str): NCBI taxonomy id. 

        Returns:
            (str): Three or four letters KEEGG organism code (Ex: hsa, lma, etc).
        """

        self.load_genome_data()

        if str(tax_id) in self.tax_ids_and_its_organisms:
            return self.tax_ids_and_its_organisms[ str(tax_id) ]
        else:
            return None



