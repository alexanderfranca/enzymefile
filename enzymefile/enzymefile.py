import re
import pprint


class EnzymeFile:
    """
    This class deals with the KEGG 'enzyme' file.

    'enzyme' file is typically stored under the ligand/ directory and represents the Enzyme database in a DBGET file format.

    'enzyme' has all entries for EC numbers, its proteins/genes associated and its metabolic pathway map numbers (among other data).

    This class is selective, in other words, not all the file data will be parsed but only that data that is more important.

    Despite 'enzyme' file has a considerable size (~150MB), in a computer with Intel Core i5, 4GB RAM, load all the file in
    memory tooks only 3:30 minutes.

    Because of that, load the entire parsed file can be more productive for searching data among the result dictionary.

    Attributes:
        enzyme_file(file): File handle of Enzyme file database.
        enzyme_entries(dict): This class is about filling this dictionary. Basically all Enzyme database data goes there.
        genes_and_its_ec_numbers(dict): Genes/protein identifications and its related EC numbers. This is an auxiliary dictionary (to make faster and simpler data querying).
        organism_and_its_ec_numbers(dict): Organism codes and its related EC numbers. This is an auxiliary dictionary (to make faster and simpler data querying).
        ec_numbers_and_its_genes(dict): EC numbers and its related gene/protein identifications. This is an auxiliary dictionary (to make faster and simpler data querying).
        genes_and_its_ec_numbers_by_organism(dict): Genes/protein identifications and its related EC numbers with organism codes as main keys. This is an auxiliary dictionary (to make faster and simpler data querying).
        genes_available(list): List of available genes/protein identifications in the Enzyme database. This is an auxiliary dictionary (to make faster and simpler data querying).
    """

    def __init__(self, reader=None):

        self.reader = reader

        # This class is about filling this list.
        self.enzyme_entries = {}

        self.genes_and_its_ec_numbers = {}
        self.organism_and_its_ec_numbers = {}
        self.ec_numbers_and_its_genes = {}
        self.genes_and_its_ec_numbers_by_organism = {}

        self.genes_available = []

    def load_enzyme_file(self):
        """
        Loads the 'enzyme' file in a dictionary into the memory, in case it's not already loaded.
        """

        if not self.is_enzyme_data_generated():
            self.generate_enzyme_data()

    def entry_entry_records(self, string=None):
        """
        Weird name for a method.

        But this method returns the data inside the entry ENTRY.

        That means, for the 'enzyme' file, the line that stores the EC number and if the EC number is obsolete or not.

        Args:
            string(str): An ENTRY string.

        Returns:
            (dict): A dictionary containing the EC number and its status.
        """

        # Separate the EC number from its followed data.
        re_ec_entry = re.compile(
            '^EC\s([0-9]{1,}\.[0-9]{1,}\.[0-9]{1,}\.[0-9]{1,}|n.*)\s(.*)$')

        result = re_ec_entry.search(string)

        # 'entry_type' is simply an string with 'Obsolete Enzyme' or simply 'Enzyme' string.
        return {'ec_number': result.group(1), 'entry_type': result.group(2)}

    def is_enzyme_data_generated(self):
        """
        Returns if the 'enzyme' file was loaded in a dictionary into the memory.
        """

        if len(self.enzyme_entries) == 0:
            return False
        else:
            return True

    def is_obsolete_ec(self, string=None):
        """
        Test if a entry ENTRY record is obsolete.

        It means an entry like 'EC 1.1.1.5        Obsolete  Enzyme'

        Args:
            string(str): The part from an entry ENTRY that has 'Obsolete Enzyme' or 'Enzyme' string.

        Returns:
            (boolean): If 'Obsolete' was found in the string.
        """

        re_obsolete_entry = re.compile('^.*Obsolete.*$')

        if re_obsolete_entry.search(string):
            return True
        else:
            return False

    def genes_entry(self, entry=None):
        """
        Returns a dictionary containing all organisms and its genes from a Enzyme entry.

        Args:
            entry(dict): Dictionary. Typically an dictionary processed by the DBGETReader class.

        Returns:
            (dict): All organisms and its genes from a Enzyme entry.
        """

        if 'GENES' in entry:
            genomes = entry['GENES']
            genome_genes = self.genome_genes(genomes)

            return genome_genes
        else:
            return {}

    def pathway_entry(self, entry=None):
        """
        Returns a list containing all metabolic pathway maps from an Enzyme entry.

        Args:
            entry(dict): Dictionary. Typically an dictionary processed by the DBGETReader class.

        Returns:
            (list): All metabolic pathways maps.
        """

        if 'PATHWAY' in entry:
            pathways = self.entry_pathways(entry['PATHWAY'])

            return pathways
        else:
            return []

    def names_entry(self, entry=None):
        """
        Returns a list containing all names for the Enzyme entry.

        Args:
            entry(dict): Dictionary. Typically an dictionary processed by the DBGETReader class.

        Returns:
            (list): All possible entry names.
        """

        if 'NAME' in entry:
            names = self.names(entry['NAME'])

            return names
        else:
            return []

    def orthology_entry(self, entry=None):
        """
        Returns a list containing all orthology ids for the Enzyme entry.

        Args:
            entry(dict): Dictionary. Typically an dictionary processed by the DBGETReader class.

        Returns:
            (list): All orthology ids.
        """

        if 'ORTHOLOGY' in entry:
            orthologies = self.orthologies(entry['ORTHOLOGY'])

            return orthologies
        else:
            return []

    # To be honest... it doesn't look like a good approach...
    def generate_genome_dictionaries(self, genome_genes=None, ec_number=None):
        """
        Generate some important class dictionaries.

        Those dictionaries make faster (more direct) to access different aspects from the genome genes entry.

        Args:
            genome_genes(dict): Dictionary. Typically a dictionary processed by the DBGETReader class.
            ec_number(str): EC Number from what this method has to extract data.

        Returns:
            (void): Only fill class dictionaries.
        """

        genome_genes = genome_genes
        ec_number = ec_number

        for genome, genes in genome_genes.items():
            for gene in genes:
                key = genome + ':' + gene

                if key not in self.genes_and_its_ec_numbers:
                    self.genes_and_its_ec_numbers[key] = []

                if genome not in self.genes_and_its_ec_numbers_by_organism:
                    self.genes_and_its_ec_numbers_by_organism[genome] = {}

                if key not in self.genes_and_its_ec_numbers_by_organism[genome]:
                    self.genes_and_its_ec_numbers_by_organism[genome][key] = []

                if ec_number not in self.ec_numbers_and_its_genes:
                    self.ec_numbers_and_its_genes[ec_number] = []

                self.genes_available.append(key)
                self.genes_and_its_ec_numbers[key].append(ec_number)
                self.genes_and_its_ec_numbers_by_organism[genome][key].append(
                    ec_number)

                self.ec_numbers_and_its_genes[ec_number].append(key)

            if genome not in self.organism_and_its_ec_numbers:
                self.organism_and_its_ec_numbers[genome] = []

            self.organism_and_its_ec_numbers[genome].append(ec_number)

    def entries_position(self):
        """
        Returns the entry positions from the file.

        Returns:
            (list): List of integers.
        """

        return self.reader.entries_position()

    def parsed_entry(self, position=None):
        """
        Returns the parsed entry in a specific file position.

        Returns:
            (dict): Entry

        """

        return self.reader.parsed_entry(position)

    def generate_enzyme_data(self):
        """
        Generate a dictionary containing enzyme data.

        This method is selective, in other words, despite the 'enzyme' file has a great amount of information,
        this method appends only some specific (and most useful) data.

        But feel free to improve the result dictionary with more fields from the 'enzyme' file.

        Returns:
            (dict): 'enzyme' file data in a dictionary format.
        """

        genome_genes = {}
        pathways = {}
        names = {}

        positions = self.entries_position()

        for position in positions:

            # Our entry. In other words, this method is all about the data
            # 'enzyme_entry' below.
            enzyme_entry = self.parsed_entry(position)

            # Make sure the entry actual exists.
            if 'ENTRY' in enzyme_entry:

                # From a DBGET file perspective, even the entry ENTRY can be more than one.
                # But we know this is the 'enzyme' KEGG file, so we know
                # there's only one.
                entry = enzyme_entry['ENTRY'][0]

                ec_entry = self.entry_entry_records(entry)

                if self.is_obsolete_ec(ec_entry['entry_type']):
                    continue

                ec_number = ec_entry['ec_number']

                # Parse the other 'enzyme' file entries.
                genome_genes = self.genes_entry(enzyme_entry)
                pathways = self.pathway_entry(enzyme_entry)
                names = self.names_entry(enzyme_entry)
                orthologies = self.orthology_entry(enzyme_entry)

                # One more dictionary.
                self.enzyme_entries[ec_number] = {
                    'ec': ec_number,
                    'gene': genome_genes,
                    'pathway': pathways,
                    'name': names,
                    'orthology': orthologies}

                # Generate some important class dictionaries.
                self.generate_genome_dictionaries(genome_genes, ec_number)

                # Reset stuff
                genome_genes = {}
                pathways = {}
                names = {}

        return self.enzyme_entries

    def orthologies(self, orthologies=None):
        """
        Returns a list of KO (orthology) numbers.

        Args:
            orthologies(list): List of string containing KO numbers and its description.

        Returns:
            (list): List of KO numbers without the descritiion.
        """

        entry_orthologies = []

        for orthology in orthologies:

            orthology_number = self.orthology_number(orthology)

            entry_orthologies.append(orthology_number)

        return entry_orthologies

    def orthology_number(self, string=None):
        """
        Returns the KO (orthology) number from a string that contains the KO number and its description.

        Args:
            string(str): The string containing KO number and its description.

        Returns:
            (str): KO number.
        """

        re_orthology = re.compile('^(K[0-9]{1,})\s(.*)$')

        result = re_orthology.search(string)
        orthology_number = result.group(1)

        return orthology_number

    def genome_genes(self, genomes=None):
        """
        Returns a dictionary containing the KEGG organism code and its protein identifications for the current EC number.

        Args:
            genomes(list): List of strings containing KEGG organism code and protein identifications.

        Returns:
            (dict): Dictionary with KEGG organism codes and its protein identifications.
        """

        genome_genes = {}

        for genome in genomes:

            organism_code = self.organism_code(genome)
            genes = self.genes(genome)

            if organism_code not in genome_genes:
                genome_genes[organism_code] = []

            genome_genes[organism_code] = genes

        return genome_genes

    def genes(self, string=None):
        """
        Returns the list of protein identification from a organism.

        Args:
            string(str): String that contains the KEGG organism code and its protein identifications.

        Returns:
            (list): List of genes (protein identifications).
        """

        result = []

        records = string.split(':')

        genes = records[1]
        genes = re.sub('^\ {1,}', '', genes)
        genes = re.sub('\ {1,}$', '', genes)
        genes = genes.split(' ')

        re_gene_code_only = re.compile('^(.*)\(.*\).*$')

        for gene in genes:

            gene = gene.lower()
            gene = gene.replace(' ', '')

            search_gene_code = re_gene_code_only.search(gene)

            if search_gene_code:
                gene = search_gene_code.group(1)

            result.append(gene)

        return result

    def entry_pathways(self, pathways=None):
        """
        Returns the metabolic pathway maps numbers.

        Args:
            pathways(list): List of strings containing map numbers and its correlated data (which will be removed).

        Returns:
            (list): List of metabolic pathway map numbers.
        """

        entry_pathways = []

        for pathway in pathways:

            map_number = self.map_number(pathway)

            entry_pathways.append(map_number)

        return entry_pathways

    def names(self, names=None):
        """
        Returns a list of names of the enzyme entry.

        Args:
            names(list): List of strings that contains the names for the enzyme.

        Returns:
            (list): List of names for the enzyme.
        """

        entry_names = []

        for name in names:

            name = self.name(name)

            entry_names.append(name)

        return entry_names

    def name(self, string=None):
        """
        Returns a string containing the name for the enzyme.

        This methos basically only treats an string to remove special chars.

        Args:
            string(str): A string line containing the name for the enzyme.

        Returns:
            (str): A name for the enzyme.
        """

        name = string.replace(';', '')
        name = re.sub('^\ {1,}', '', name)
        name = re.sub('\ {1,}$', '', name)

        return name

    def map_number(self, string=None):
        """
        Returns the metabolic pathway map number from a string.

        Args:
            string(str): An string containing a map number.

        Returns:
            (str): The map number without any useless char/string.
        """

        re_pathway = re.compile('^ec([0-9]{1,})\s(.*)$')

        result = re_pathway.search(string)
        map_number = result.group(1)

        return map_number

    def organism_code(self, string=None):
        """
        Returns the organism code from a string.

        Args:
            string(str): An string containing the organism code.

        Returns:
            (str): The organism code in lower case and withou any special or useless char.
        """

        records = string.split(':')

        organism_code = records[0]
        organism_code = organism_code.lower()
        organism_code = organism_code.replace(' ', '')

        return organism_code

    def maps_by_ec_number(self, ec_number=None):
        """
        Returns all the metabolic pathway map number from a EC number.

        Args:
            ec_number(str): The EC number.

        Returns:
            (list): List of metabolic pathway maps number.
        """

        self.load_enzyme_file()

        if ec_number in self.enzyme_entries:
            return self.enzyme_entries[ec_number]['pathway']
        else:
            return []

    def names_by_ec_number(self, ec_number=None):
        """
        Returns all the names from a EC number.

        Args:
            ec_number(str): The EC number.

        Returns:
            (list): List of names.
        """

        self.load_enzyme_file()

        return self.enzyme_entries[ec_number]['name']

    def organisms_by_ec_number(self, ec_number=None):
        """
        Returns all KEGG organisms code from a EC number.

        Args:
            ec_number(str): The EC number.

        Returns:
            (list): List of KEGG organism codes.
        """

        self.load_enzyme_file()

        organisms = []

        genes = self.enzyme_entries[ec_number]['gene']

        for g in genes:
            organisms.append(g)

        return organisms

    def ec_numbers_by_organism(self, organism_code=None):
        """
        Returns all EC numbers associated with the three/four letter KEGG orgnanism code (Ex: hsa, lma, etc)

        Args:
            organism_code(str): Three/four letter KEGG organism code.

        Returns:
            (list): EC numbers.
        """

        self.load_enzyme_file()

        if organism_code in self.organism_and_its_ec_numbers:
            return self.organism_and_its_ec_numbers[organism_code]
        else:
            return []

    def organism_and_genes_by_ec_number(self, ec_number=None):
        """
        Returns a dictionary containing KEGG organism codes and its protein identifications by EC number.

        Args:
            ec_number(str): EC number.

        Returns:
            (dict): Organisms and protein identifications for the EC number.
        """

        organism_and_genes = {}

        self.load_enzyme_file()

        genes = self.enzyme_entries[ec_number]['gene']

        for organism, ids in genes.items():

            organism_and_genes[organism] = ids

        return organism_and_genes

    def ec_numbers_by_gene(self, gene=None):
        """
        Returns all the EC numbers of a gene.

        Args:
            gene(str): Gene identification in a format: organism code + : + identification. Ex: 'xla:108713048'.

        Returns:
            (list): EC numbers.
        """

        self.load_enzyme_file()

        ecs = []

        if gene in self.genes_available:
            ecs = self.genes_and_its_ec_numbers[gene]

        return ecs

    def ec_numbers_by_gene_and_organism(self, gene=None, organism=None):
        """
        Return EC numbers from a gene/protein identification and specific organism.

        Args:
            gene(str): Gene/protein identification. Example: 'hsa:115201'.
            organism(str): KEGG 3 or 4 letters organism code. Example: 'hsa'.

        Returns:
            (list): EC numbers.
        """

        self.load_enzyme_file()

        data = {}

        if organism in self.genes_and_its_ec_numbers_by_organism:
            if gene in self.genes_and_its_ec_numbers_by_organism[organism]:
                data = self.genes_and_its_ec_numbers_by_organism[organism][gene]
        else:
            data = {}

        return data

    def all_genes_and_its_ec_numbers(self):
        """
        Returns a dictionary containing all the genes and its related EC numbers.

        Returns:
            (dict): Genes and EC numbers.
        """

        self.load_enzyme_file()

        return self.genes_and_its_ec_numbers

    def genes_by_ec_number(self, ec_number=None):
        """
        Returns all gene identifications associated with the EC number.

        Args:
            ec_number(str): EC number.

        Returns:
            (list): Gene identifications.
        """

        self.load_enzyme_file()

        if ec_number in self.ec_numbers_and_its_genes:
            return self.ec_numbers_and_its_genes[ec_number]
        else:
            return []

    def complete_ec_numbers(self):
        """
        Returns the list of all complete EC numbers (ECs with all 4 digitis fully filled).

        Returns:
            (list): EC numbers
        """

        self.load_enzyme_file()

        complete_ecs = []

        re_incomplete = re.compile('-')

        for ec in self.enzyme_entries.keys():
            result = re_incomplete.search(ec)

            if not result:
                complete_ecs.append(ec)

        complete_ecs = set(complete_ecs)
        complete_ecs = list(complete_ecs)

        return complete_ecs

    def incomplete_ec_numbers(self):
        """
        Returns the list of all INcomplete EC numbers (ECs with some of its 4 digitis filled with dash '-').

        Returns:
            (list): EC numbers
        """

        self.load_enzyme_file()

        incomplete_ecs = []

        re_incomplete = re.compile('-')

        for ec in self.enzyme_entries.keys():
            result = re_incomplete.search(ec)

            if result:
                incomplete_ecs.append(ec)

        incomplete_ecs = set(incomplete_ecs)
        incomplete_ecs = list(incomplete_ecs)

        return incomplete_ecs

    def all_ec_numbers(self):
        """
        Returns a list containing all the EC numbers from the Enzyme file.
        """

        self.load_enzyme_file()

        complete = self.complete_ec_numbers()
        incomplete = self.incomplete_ec_numbers()

        all_ecs = complete + incomplete

        return list(all_ecs)

    def all_organism_ecs(self):
        """
        Return all organism codes and its related EC numbers.

        Returns:
            (dict): KEGG organism codes as keys and EC numbers as list of values.
        """

        self.load_enzyme_file()

        return self.organism_and_its_ec_numbers
