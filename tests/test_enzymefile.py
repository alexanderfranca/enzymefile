import sys
import os
import unittest
from enzymefile.enzymefile import EnzymeFile 
from dbgetreader import *
import re
import pprint

class test_EnzymeKegg(unittest.TestCase):

    def setUp(self):

        dbget = DBGET(file_to_read='./tests/fixtures/enzyme')
        dbgetr = DBGETReader(reader=dbget)

        self.enzyme = EnzymeFile(dbgetr)

    def test_entry_entry_records(self):

        string = 'EC 1.1.1.5        Obsolete  Enzyme'

        data = self.enzyme.entry_entry_records(string)

        self.assertEqual(data['ec_number'], '1.1.1.5')

    def test_is_obsolete_ec(self):

        string = 'EC 1.1.1.5        Obsolete  Enzyme'

        data = self.enzyme.entry_entry_records(string)

        self.assertTrue(self.enzyme.is_obsolete_ec(data['entry_type']))

    def test_complete_ec_numbers(self):

        result = self.enzyme.complete_ec_numbers()

        expected = [ '1.1.1.1', '2.7.7.96', '3.6.1.13', '1.1.1.2', '1.1.1.3', '3.6.1.59', '3.6.1.58' ]

        self.assertEqual( len(result), 7 )
        self.assertTrue( set(expected) == set(result) )

    def test_incomplete_ec_numbers(self):

        result = self.enzyme.incomplete_ec_numbers()

        expected = []

        self.assertEqual( len(result), 0 )
        self.assertTrue( set(expected) == set(result) )

    def test_generate_enzyme_data(self):

        self.enzyme.generate_enzyme_data()

        for enzyme in self.enzyme.enzyme_entries:
            data = self.enzyme.enzyme_entries[enzyme]
            break

        read_keys = data.keys()

        expected_keys = ['ec', 'gene', 'pathway', 'name', 'orthology']

        if len(
                set(expected_keys) -
                set(read_keys)) > 0 or len(
                set(read_keys) -
                set(expected_keys)) > 0:
            keys_changed = True
        else:
            keys_changed = False

        self.assertFalse(keys_changed)

    def test_orthologies(self):

        self.enzyme.generate_enzyme_data()

        for enzyme in self.enzyme.enzyme_entries:
            data = self.enzyme.enzyme_entries[enzyme]
            break

        self.assertEqual(data['orthology'][0], 'K00001')

    def test_orthology_number(self):

        string = 'K07514  enoyl-co_a hydratase / 3-hydroxyacyl-co_a dehydrogenase / 3,2-trans-enoyl-co_a isomerase'

        self.assertEqual(self.enzyme.orthology_number(string), 'K07514')

    def test_genome_genes(self):

        self.enzyme.entries_position()
        data = self.enzyme.parsed_entry(0)

        names = self.enzyme.names(data['NAME'])

        self.assertTrue(isinstance(names, list))
        self.assertTrue(len(names) > 0)

    def test_entry_pathways(self):

        self.enzyme.entries_position()
        data = self.enzyme.parsed_entry(0)

        pathways = self.enzyme.entry_pathways(data['PATHWAY'])

        self.assertTrue(isinstance(pathways, list))
        self.assertTrue(len(pathways) > 0)

    def test_name(self):

        string = '    NAD-specific aromatic alcohol dehydrogenase;  '

        self.assertEqual(self.enzyme.name(string),
                          'NAD-specific aromatic alcohol dehydrogenase')

    def test_map_number(self):

        string = 'ec00260  Glycine, serine and threonine metabolism'

        self.assertEqual(self.enzyme.map_number(string), '00260')

    def test_organism_code(self):

        string = 'HSA: 124(ADH1A) 125(ADH1B) 126(ADH1C) 127(ADH4) 128(ADH5) 130(ADH6) 131(ADH7)'

        self.assertEqual(self.enzyme.organism_code(string), 'hsa')

    def test_maps_by_ec_number(self):

        expected_maps = [
            '00010',
            '00071',
            '00260',
            '00350',
            '00592',
            '00625',
            '00626',
            '00830',
            '00980',
            '00982',
            '01100',
            '01110',
            '01120',
            '01130']

        self.assertEqual(
            self.enzyme.maps_by_ec_number('1.1.1.1'),
            expected_maps)

    def test_names_by_ec_number(self):

        expected_names = ['homoserine dehydrogenase', 'HSDH', 'HSD']

        self.assertEqual(
            self.enzyme.names_by_ec_number('1.1.1.3'),
            expected_names)

    def test_organisms_by_ec_number(self):

        self.max_diff = None

        self.assertEqual(
            len(self.enzyme.organisms_by_ec_number('1.1.1.1')), 3427)

    def test_organisms_and_genes_by_ec_number(self):

        self.max_diff = None

        self.assertEqual(
            len(self.enzyme.organism_and_genes_by_ec_number('1.1.1.1')), 3427)

    def test_ec_numbers_by_organism(self):

        ec_numbers = self.enzyme.ec_numbers_by_organism('hsa')

        expected_ecs = [
            '1.1.1.1',
            '1.1.1.2',
            '1.1.1.3',
            '3.6.1.13',
            '3.6.1.58',
            '3.6.1.59',
            '2.7.7.96']

        self.assertEqual(ec_numbers, expected_ecs)

    def test_ec_numbers_by_organism_void(self):

        ec_numbers = self.enzyme.ec_numbers_by_organism(
            'this_organism_doesnt_exist')
        self.assertEqual(ec_numbers, [])

    def test_organism_and_genes_by_ec_number(self):

        data = self.enzyme.organism_and_genes_by_ec_number('2.7.7.96')

        self.assertTrue(len(data) == 9)
        self.assertEqual(data['xla'], ['108713048', '432263'])

    def test_ec_numbers_by_gene(self):

        data = self.enzyme.ec_numbers_by_gene('xla:108713048')

        self.assertEqual(data, ['3.6.1.13', '3.6.1.58', '2.7.7.96'])

    def test_ec_numbers_by_gene_void(self):

        data = self.enzyme.ec_numbers_by_gene('theres_not_found_gene_here')

        self.assertEqual(data, [])

    def test_all_genes_and_its_ec_numbers(self):

        data = self.enzyme.all_genes_and_its_ec_numbers()

        self.assertTrue(isinstance(data, dict))
        self.assertEqual(
            data['xla:108713048'], [
                '3.6.1.13', '3.6.1.58', '2.7.7.96'])

    def test_genes_by_ec_number(self):

        data = self.enzyme.genes_by_ec_number('2.7.7.96')

        expected_genes = [
            'hsa:11164',
            'ptr:450306',
            'pps:100977876',
            'ggo:101148720',
            'xla:108713048',
            'xla:432263',
            'xtr:549367',
            'sasa:100196024',
            'sasa:106609278',
            'acan:aca1_131000',
            'acan:aca1_177470',
            'gtt:guithdraft_166567']

        self.assertEqual(data, expected_genes)

    def test_genes_by_ec_number_void(self):

        data = self.enzyme.genes_by_ec_number('this_ec_number_doesnt_exist')

        expected_genes = []

        self.assertEqual(data, expected_genes)

    def test_all_ec_numbers(self):

        data = self.enzyme.all_ec_numbers()

        self.assertTrue(len(data) > 1)
        self.assertTrue(isinstance(data, list))

    def test_genes_entry(self):

        positions = self.enzyme.entries_position()

        for position in positions:
            entry = self.enzyme.parsed_entry(position)
            break

        genes_entry = self.enzyme.genes_entry(entry)

        self.assertTrue(isinstance(genes_entry, dict))
        self.assertTrue(len(genes_entry) > 1)

    def test_pathway_entry(self):

        positions = self.enzyme.entries_position()

        for position in positions:
            entry = self.enzyme.parsed_entry(position)
            break

        data = self.enzyme.pathway_entry(entry)

        self.assertTrue(isinstance(data, list))
        self.assertTrue(len(data) > 1)

    def test_names_entry(self):

        positions = self.enzyme.entries_position()

        for position in positions:
            entry = self.enzyme.parsed_entry(position)
            break

        data = self.enzyme.names_entry(entry)

        self.assertTrue(isinstance(data, list))
        self.assertTrue(len(data) > 1)

    def test_orthology_entry(self):

        positions = self.enzyme.entries_position()

        for position in positions:
            entry = self.enzyme.parsed_entry(position)
            break

        data = self.enzyme.orthology_entry(entry)

        self.assertTrue(isinstance(data, list))
        self.assertTrue(len(data) > 1)

    def test_ec_numbers_by_gene_and_organism(self):

        gene = 'hsa:115201'
        organism = 'hsa'

        data = self.enzyme.ec_numbers_by_gene_and_organism(gene, organism)

        expected_ecs = ['1.1.1.1', '3.6.1.59']

        self.assertTrue(len(data) > 1)
        self.assertTrue(isinstance(data, list))
        self.assertEqual(data, expected_ecs)

    def test_all_organism_ecs(self):

        data = self.enzyme.all_organism_ecs()

        expected_ecs = [
            '1.1.1.1',
            '1.1.1.2',
            '1.1.1.3',
            '3.6.1.13',
            '3.6.1.59']

        self.assertTrue(len(data) > 1)
        self.assertTrue(isinstance(data, dict))
        self.assertEqual(data['abe'], expected_ecs)


if __name__ == "__main__":
    unittest.main()
