import re
import argparse
import os
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
import Bio
from collections import OrderedDict
import copy


# 找到所有的基因
interval = 0

def find_trans(records: [SeqIO.SeqRecord]):
    trans = []
    for record in records:
        for feature in record.features:
            if feature.type == 'CDS':
                if 'exception' in feature.qualifiers.keys():
                    if feature.qualifiers['exception'][0] == 'trans-splicing':
                        gene_name = feature.qualifiers['gene'][0]
                        while gene_name in trans:
                            gene_name += ' '
                        trans.append(gene_name)
    return trans

def find_genes(record: SeqIO.SeqRecord):
    feature_types_to_shorten = ['CDS', 'gene', 'tRNA', 'rRNA', 'intron']
    genes = []  # gene位置
    # 找到feature_types_to_shorten中类型的feature
    for feature in record.features:
        if feature.type in feature_types_to_shorten:
            for part in feature.location.parts:
                if (part.nofuzzy_start, part.nofuzzy_end) not in genes:
                    genes.append((part.nofuzzy_start, part.nofuzzy_end))
    return genes


def move_gene(location, genes, position_bias):
    if isinstance(location, FeatureLocation):
        index = genes.index((location.nofuzzy_start, location.nofuzzy_end))
        return FeatureLocation(genes[index][0] - position_bias[index],
                               genes[index][1] - position_bias[index],
                               ref=location.ref,
                               ref_db=location.ref_db,
                               strand=location.strand)
    elif isinstance(location, CompoundLocation):
        for j in range(len(location.parts)):
            location.parts[j] = move_gene(location.parts[j], genes, position_bias)
        return location
    else:
        raise


def decouple_location(location):
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


exon_type = ['tRNA', 'rRNA', 'other']


class genes_dict(dict):
    def __init__(self):
        super(genes_dict, self).__init__()

    def add(self, gene_name, feature_name, feature):
        if gene_name not in self.keys():
            self[gene_name] = {}
        if feature_name == 'exon':
            if 'exon' not in self[gene_name].keys():
                self[gene_name]['exon'] = {}
            if feature.qualifiers['number'][0] in self[gene_name]['exon'].keys():
                self.add(gene_name + ' ', feature_name, feature)
            else:
                self[gene_name]['exon'][feature.qualifiers['number'][0]] = feature
        else:
            if feature_name in self[gene_name].keys():
                self.add(gene_name + ' ', feature_name, feature)
            else:
                self[gene_name][feature_name] = feature


def shorten_distance(records: [SeqIO.SeqRecord], ratio, out=''):
    trans = find_trans(records)
    # 平移features
    source_feature = records[0].features[0]
    genes = genes_dict()
    for i in range(len(records)):
        for feature in records[i].features:
            if feature.type in ['CDS', 'genes', 'exon']:
                gene_name = feature.qualifiers['gene'][0]
                genes.add(gene_name, feature.type, feature)

    for gene_name in genes.keys():
        if 'exon' in genes[gene_name].keys():
            if len(genes[gene_name]['exon']) > 1:
                start = interval
                new_features = [source_feature]

                for i in range(len(genes[gene_name]['exon'].values())):
                    for feature in genes[gene_name]['exon'].values():
                        j = int(feature.qualifiers['number'][0])
                        if j == i+1:
                            location = feature.location
                            # feature = copy.deepcopy(genes[gene_name]['CDS'])
                            new_location = FeatureLocation(start,
                                                           start + 1 + location.nofuzzy_end - location.nofuzzy_start,
                                                           ref=location.ref,
                                                           ref_db=location.ref_db,
                                                           strand=location.strand)
                            new_location.strand = 1

                            new_features.append(SeqFeature(location=new_location, type='gene',
                                                           qualifiers=OrderedDict({'gene': [f'exon{j}({location.nofuzzy_end - location.nofuzzy_start})']})))
                            feature.qualifiers['gene'] = [f'exon{j}({location.nofuzzy_end - location.nofuzzy_start})']
                            feature.location = new_location
                            feature.type = 'exon' + feature.qualifiers['number'][0]  # exon_type[j % len(exon_type)]
                            new_features.append(feature)
                            start = start + 1 + location.nofuzzy_end - location.nofuzzy_start + interval
                            break
                        # j += 1
                new_record = copy.deepcopy(records[0])
                new_record.seq = records[0].seq[:start]
                new_features[0].location = FeatureLocation(0,
                                                           start,
                                                           ref=new_features[0].location.ref,
                                                           ref_db=new_features[0].location.ref_db,
                                                           strand=new_features[0].location.strand)
                new_record.features = new_features
                new_record = shorten_distance1(new_record, ratio)
                if gene_name in trans:
                    with open(os.path.join(out, 'trans', gene_name + '.gb'), 'w') as f:
                        myGenBankWriter(f).write_file([new_record])
                else:
                    with open(os.path.join(out, 'cis', gene_name + '.gb'), 'w') as f:
                        myGenBankWriter(f).write_file([new_record])

    # while i < len(records.features):
    #     if records.features[i].type == 'CDS':
    #         locations = decouple_location(records.features[i].location)
    #         if len(locations) > 1:
    #             start = 200
    #             new_features = [source_feature]
    #             j = 1
    #             gene_name = records.features[i].qualifiers['gene'][0]
    #             for location in locations:
    #                 feature = copy.deepcopy(records.features[i])
    #                 new_location = FeatureLocation(start,
    #                                                start + 1 + location.nofuzzy_end - location.nofuzzy_start,
    #                                                ref=location.ref,
    #                                                ref_db=location.ref_db,
    #                                                strand=location.strand)
    #                 new_location.strand = 1
    #
    #                 new_features.append(SeqFeature(location=new_location, type='gene',
    #                                                qualifiers=OrderedDict({'gene': [f'exon{j}({location.nofuzzy_end - location.nofuzzy_start})']})))
    #                 feature.qualifiers['gene'] = [f'exon{j}({location.nofuzzy_end - location.nofuzzy_start})']
    #                 feature.location = new_location
    #                 feature.type = exon_type[j % len(exon_type)]
    #                 new_features.append(feature)
    #                 start = start + 1 + location.nofuzzy_end - location.nofuzzy_start + 200
    #                 j += 1
    #             new_record = copy.deepcopy(records)
    #             new_record.seq = records.seq[:start]
    #             new_features[0].location = FeatureLocation(0,
    #                                                        start,
    #                                                        ref=new_features[0].location.ref,
    #                                                        ref_db=new_features[0].location.ref_db,
    #                                                        strand=new_features[0].location.strand)
    #             new_record.features = new_features
    #             new_record = shorten_distance1(new_record, ratio)
    #             if gene_name in ['nad1', 'nad2', 'nad5', 'cox2', 'rps3', 'nad4']:
    #                 with open(os.path.join(out, 'trans', gene_name + '.gb'), 'w') as f:
    #                     myGenBankWriter(f).write_file([new_record])
    #             else:
    #                 with open(os.path.join(out, 'cis', gene_name + '.gb'), 'w') as f:
    #                     myGenBankWriter(f).write_file([new_record])
    #         i += 1
    #     else:
    #         i += 1
    return


def shorten_distance1(record: SeqIO.SeqRecord, ratio):
    genes = find_genes(record)
    genes.sort()
    bias_ratio = 1 - ratio
    original_length = len(record.seq)

    # 计算向前平移的距离
    position_bias = []
    last_end = 0
    bias = 0
    for i in range(len(genes)):
        if genes[i][0] > last_end:
            bias += int((genes[i][0] - last_end) * bias_ratio)
        last_end = max(genes[i][1], last_end)
        position_bias.append(bias)

    # 找出初始文件中基因位置的最大值，算出最大值与总长度的比例，用该比例与修改后基因位置的最大值得出修改后总长度
    if len(position_bias) > 0:
        new_length = int((original_length - last_end) * ratio + last_end - bias)
    else:
        new_length = original_length

    # 平移features
    i = 0
    new_features = []
    while i < len(record.features):
        if record.features[i].type == 'source':
            feature = copy.deepcopy(record.features[i])  # 因为后面这个record还要用到，我们要避免对record的直接修改
            feature.location = FeatureLocation(0, new_length, strand=feature.location.strand)
            new_features.append(feature)
            i += 1
        else:
            location = record.features[i].location
            feature = copy.deepcopy(record.features[i])
            new_location = move_gene(location, genes, position_bias)  # 缩进location
            feature.location = new_location
            new_features.append(feature)
            i += 1
    record.features = new_features
    record.seq = record.seq[:new_length]
    return record

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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='the RNA file to edit', type=str, nargs='+', default=['.\Saccharum_officinarum_chr1.gb','.\Saccharum_officinarum_chr2.gb'])
    parser.add_argument('--out', help='the path to save result', type=str, default='dna3')
    parser.add_argument('--ratio', help='the ratio of interval to be kept', type=float, default=0.3)
    args = parser.parse_args()
    if not os.path.isdir(os.path.join(args.out, 'trans')):
        os.makedirs(os.path.join(args.out, 'trans'))
    if not os.path.isdir(os.path.join(args.out, 'cis')):
        os.makedirs(os.path.join(args.out, 'cis'))
    records = [SeqIO.read(file_name, "genbank") for file_name in args.input]  # 读取gb文件
    for i in range(len(records)):
        records[i].annotations['organism'] = records[i].annotations['organism'].replace(' ', '\n            ', 1)
    shorten_distance(records, args.ratio, args.out)
