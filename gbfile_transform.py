import argparse
import copy
import os
from collections import OrderedDict
import Bio
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from my_error import *

def decouple_location(location):
    """
    把location分解成单个外显子组成的列表
    """
    if isinstance(location, FeatureLocation):
        return [location]
    elif isinstance(location, CompoundLocation):
        locations = []
        for i in location.parts:
            locations += decouple_location(i)
        return locations
    elif isinstance(location, list):
        locations = []
        for i in location:
            locations += decouple_location(i)
        return locations
    else:
        raise


from Bio import BiopythonWarning
import warnings


class myGenBankWriter(SeqIO.InsdcIO.GenBankWriter):
    def _write_single_line(self, tag, text, flag=True):
        """Write single line in each GenBank record (PRIVATE).

        Used in the 'header' of each GenBank record.
        """
        assert len(tag) < self.HEADER_WIDTH
        if len(text) > self.MAX_WIDTH - self.HEADER_WIDTH:
            if tag:
                warnings.warn(
                    "Annotation %r too long for %r line" % (text, tag), BiopythonWarning
                )
            else:
                # Can't give such a precise warning
                warnings.warn("Annotation %r too long" % text, BiopythonWarning)
        if flag:
            self.handle.write(
                "%s%s\n" % (tag.ljust(self.HEADER_WIDTH), text.replace("\n", " "))
            )
        else:
            self.handle.write(
                "%s%s\n" % (tag.ljust(self.HEADER_WIDTH), text)
            )

    def write_record(self, record):
        """Write a single record to the output file."""
        handle = self.handle
        self._write_the_first_line(record)

        default = record.id
        if default.count(".") == 1 and default[default.index(".") + 1:].isdigit():
            # Good, looks like accession.version and not something
            # else like identifier.start-end
            default = record.id.split(".", 1)[0]
        accession = self._get_annotation_str(
            record, "accession", default, just_first=True
        )
        acc_with_version = accession
        if record.id.startswith(accession + "."):
            try:
                acc_with_version = "%s.%i" % (
                    accession,
                    int(record.id.split(".", 1)[1]),
                )
            except ValueError:
                pass
        gi = self._get_annotation_str(record, "gi", just_first=True)

        descr = record.description
        if descr == "<unknown description>":
            descr = ""  # Trailing dot will be added later

        # The DEFINITION field must end with a period
        # see ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt [3.4.5]
        # and discussion https://github.com/biopython/biopython/pull/616
        # So let's add a period
        descr += "."
        self._write_multi_line("DEFINITION", descr)

        self._write_single_line("ACCESSION", accession)
        if gi != ".":
            self._write_single_line("VERSION", "%s  GI:%s" % (acc_with_version, gi))
        else:
            self._write_single_line("VERSION", "%s" % acc_with_version)

        # The NCBI initially expected two types of link,
        # e.g. "Project:28471" and "Trace Assembly Archive:123456"
        #
        # This changed and at some point the formatting switched to
        # include a space after the colon, e.g.
        #
        # LOCUS       NC_000011               1606 bp    DNA     linear   CON 06-JUN-2016
        # DEFINITION  Homo sapiens chromosome 11, GRCh38.p7 Primary Assembly.
        # ACCESSION   NC_000011 REGION: complement(5225466..5227071) GPC_000001303
        # VERSION     NC_000011.10  GI:568815587
        # DBLINK      BioProject: PRJNA168
        #             Assembly: GCF_000001405.33
        # ...
        #
        # Or,
        #
        # LOCUS       JU120277                1044 bp    mRNA    linear   TSA 27-NOV-2012
        # DEFINITION  TSA: Tupaia chinensis tbc000002.Tuchadli mRNA sequence.
        # ACCESSION   JU120277
        # VERSION     JU120277.1  GI:379775257
        # DBLINK      BioProject: PRJNA87013
        #             Sequence Read Archive: SRR433859
        # ...
        dbxrefs_with_space = []
        for x in record.dbxrefs:
            if ": " not in x:
                x = x.replace(":", ": ")
            dbxrefs_with_space.append(x)
        self._write_multi_entries("DBLINK", dbxrefs_with_space)
        del dbxrefs_with_space

        try:
            # List of strings
            # Keywords should be given separated with semi colons,
            keywords = "; ".join(record.annotations["keywords"])
            # with a trailing period:
            if not keywords.endswith("."):
                keywords += "."
        except KeyError:
            # If no keywords, there should be just a period:
            keywords = "."
        self._write_multi_line("KEYWORDS", keywords)

        if "segment" in record.annotations:
            # Deal with SEGMENT line found only in segmented records,
            # e.g. AH000819
            segment = record.annotations["segment"]
            if isinstance(segment, list):
                assert len(segment) == 1, segment
                segment = segment[0]
            self._write_single_line("SEGMENT", segment)

        self._write_multi_line("SOURCE", self._get_annotation_str(record, "source"))
        # The ORGANISM line MUST be a single line, as any continuation is the taxonomy
        org = self._get_annotation_str(record, "organism")
        if len(org) > self.MAX_WIDTH - self.HEADER_WIDTH:
            org = org[: self.MAX_WIDTH - self.HEADER_WIDTH - 4] + "..."
        self._write_single_line("  ORGANISM", org, False)
        try:
            # List of strings
            # Taxonomy should be given separated with semi colons,
            taxonomy = "; ".join(record.annotations["taxonomy"])
            # with a trailing period:
            if not taxonomy.endswith(".") and taxonomy != '':
                taxonomy += "."
        except KeyError:
            taxonomy = "."
        if taxonomy != '':
            self._write_multi_line("", taxonomy)

        if "db_source" in record.annotations:
            # Hack around the issue of BioSQL loading a list for the db_source
            db_source = record.annotations["db_source"]
            if isinstance(db_source, list):
                db_source = db_source[0]
            self._write_single_line("DBSOURCE", db_source)

        if "references" in record.annotations:
            self._write_references(record)

        if (
                "comment" in record.annotations
                or "structured_comment" in record.annotations
        ):
            self._write_comment(record)

        handle.write("FEATURES             Location/Qualifiers\n")
        rec_length = len(record)
        for feature in record.features:
            self._write_feature(feature, rec_length)
        self._write_sequence(record)
        handle.write("//\n")


def transform_gb(input_file, output):
    # print(input_file)
    records = [SeqIO.read(file_name, "genbank") for file_name in input_file]  # 读取gb文件
    genes = {}
    type1 = ['CDS', 'tRNA', 'rRNA']  # 需要提取exon的feature类型
    for j in range(len(records)):
        gene_names = {'gene': [], 'product': [], 'locus_tag': []}
        for i in range(1, len(records[j].features)):
            # if records[j].features[i].type in type1 or records[j].features[i].type == 'gene' or records[j].features[
                # i].type == 'exon':
            if records[j].features[i].type in type1 or records[j].features[i].type == 'gene':
                # assert 'gene' in records[j].features[i].qualifiers.keys() \
                #        or 'product' in records[j].features[i].qualifiers.keys() \
                #        or 'locus_tag' in records[j].features[i].qualifiers.keys(), \
                #     f'there are no gene, product, or locus_tag in feature \n{records[j].features[i]}'
                if not ('gene' in records[j].features[i].qualifiers.keys()
                        or 'product' in records[j].features[i].qualifiers.keys()
                        or 'locus_tag' in records[j].features[i].qualifiers.keys()):
                    raise noGeneNameError(records[j].features[i])
                index = -1
                for key in ['gene', 'product', 'locus_tag']:
                    if key in records[j].features[i].qualifiers.keys():
                        if records[j].features[i].qualifiers[key][0] in gene_names[key]:
                            index = gene_names[key].index(records[j].features[i].qualifiers[key][0])
                            for key2 in ['gene', 'product', 'locus_tag']:
                                if key2 in records[j].features[i].qualifiers.keys():
                                    if gene_names[key2][index] is not None and gene_names[key2][index] != records[j].features[i].qualifiers[key2][0]:
                                        index = -1
                                        break
                            if index != -1:
                                break
                if index == -1:
                    for key in ['gene', 'product', 'locus_tag']:
                        gene_names[key].append(None)
                for key in ['gene', 'product', 'locus_tag']:
                    if key in records[j].features[i].qualifiers.keys():
                        if gene_names[key][index] is None:
                            gene_names[key][index] = records[j].features[i].qualifiers[key][0]
                        elif gene_names[key][index] != records[j].features[i].qualifiers[key][0]:
                            for key2 in ['gene', 'product', 'locus_tag']:
                                gene_names[key2].append(gene_names[key2][index])
                            index = -1
                            gene_names[key][index] = records[j].features[i].qualifiers[key][0]
        for i in range(1, len(records[j].features)):
            # if (records[j].features[i].type in type1 or records[j].features[i].type == 'gene'
            #     or records[j].features[i].type == 'exon') and 'gene' not in records[j].features[i].qualifiers.keys():
            if (records[j].features[i].type in type1 or records[j].features[i].type == 'gene') and 'gene' not in records[j].features[i].qualifiers.keys():
                index = -1
                for key in ['product', 'locus_tag']:
                    if key in records[j].features[i].qualifiers.keys():
                        index = gene_names[key].index(records[j].features[i].qualifiers[key][0])
                        break
                assert index != -1, 'An unexpected error occurred in the program: information about gene name is not fully collected.'
                for key in ['gene', 'product', 'locus_tag']:
                    if gene_names[key][index] is not None:
                        records[j].features[i].qualifiers['gene'] = [gene_names[key][index]]
                        break
            if type(records[j].features[i].location) == Bio.SeqFeature.CompoundLocation:
                for h in range(len(records[j].features[i].location.parts) - 1, -1, -1):
                    if records[j].features[i].location.parts[h].ref is not None:
                        records[j].features[i].location.parts.pop(h)
            if records[j].features[i].type in type1:
                gene_name = records[j].features[i].qualifiers['gene'][0]
                while gene_name in genes.keys() and (
                        records[j].features[i].type != 'CDS' or j in genes[gene_name]['file']):
                    gene_name += ' '
                if gene_name not in genes.keys():
                    genes[gene_name] = {'location': [], 'file': [], 'cis_flag': True, 'exon_flag': False, 'type': None}
                if genes[gene_name]['type'] == None and records[j].features[i].type in type1:
                    genes[gene_name]['type'] = records[j].features[i].type
                genes[gene_name]['file'].append(j)
                locations = decouple_location(records[j].features[i].location)
                for location in locations:
                    if location not in genes[gene_name]['location']:
                        genes[gene_name]['location'].append(location)
                # if j not in genes[gene_name]['file']:
                #     genes[gene_name]['file'].append(j)
            # elif records[j].features[i].type == 'exon':
            #     gene_name = records[j].features[i].qualifiers['gene'][0]
            #     if gene_name not in genes.keys():
            #         genes[gene_name] = {'location': [], 'file': [], 'cis_flag': True, 'exon_flag': False, 'type': None}
            #     genes[gene_name]['exon_flag'] = True  # 已标注外显子就不再对该gene做处理

    # 判断cis还是trans
    for i in genes.keys():
        if genes[i]['type'] == 'CDS':
            if len(genes[i]['file']) > 1:  # 分布在多个文件上就是trans
                genes[i]['cis_flag'] = False
            else:
                assert len(genes[i]['location']) > 0, f'gene {i} only has gene or exon!!!'
                location_flag = genes[i]['location'][0].strand
                for j in range(len(genes[i]['location']) - 1):
                    if genes[i]['location'][j].strand != location_flag:  # 正反链均有外显子分布就是trans
                        genes[i]['cis_flag'] = False
                    elif (int(genes[i]['location'][j].end) - int(genes[i]['location'][j + 1].start)) * genes[i]['location'][j].strand > 0:
                            genes[i]['cis_flag'] = False

    for j in range(len(records)):
        new_features = []
        genes_in_j = []
        for i in range(len(records[j].features)):
            if records[j].features[i].type in type1:
                gene_name = records[j].features[i].qualifiers['gene'][0]
                while gene_name in genes_in_j:
                    gene_name += ' '
                genes_in_j.append(gene_name)
                if genes[gene_name]['cis_flag']:
                    if len(genes[gene_name]['location']) > 1:
                        records[j].features[i].qualifiers['exception'] = ['cis-splicing']
                else:
                    records[j].features[i].qualifiers['exception'] = ['trans-splicing']
                new_features.append(records[j].features[i])
                if len(genes[gene_name]['location']) >= 2 and genes[gene_name]['exon_flag'] == False:
                    locations = decouple_location(records[j].features[i].location)
                    for location in locations:
                        new_features.append(SeqFeature(location=location, type='exon',
                                                       qualifiers=OrderedDict({'gene': [gene_name.replace(' ', '')],
                                                                               'number': [
                                                                                   genes[gene_name]['location'].index(
                                                                                       location) + 1]})))
            # else:
            elif records[j].features[i].type != 'exon':
                new_features.append(records[j].features[i])
        new_record = copy.deepcopy(records[j])
        new_record.features = new_features
        with open(output[j], 'w') as f:
            myGenBankWriter(f).write_file([new_record])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='the RNA file to edit', type=str, default='NC_002693.gb')
    parser.add_argument('--output', help='the RNA file to edit', type=str, default='test.gb')
    # parser.add_argument('--input', help='the RNA file to edit', type=str, nargs='+',
    #                    default=['.\coffee1.gb', '.\coffee2.gb'])
    # parser.add_argument('--output', help='the RNA file to edit', type=str, nargs='+',
    #                    default=['.\coffee1_1.gb', '.\coffee2_1.gb'])
    args = parser.parse_args()
    args.input = args.input.split(',')
    args.output = args.output.split(',')
    transform_gb(args.input, args.output)
