import sys
import os

import circle_draw_ch_new
import linear_draw_ch_new

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import circle_draw_new, linear_draw_new, dna2, dna3, new_drawgenemap, \
    my_cis_and_trans, my_cis_and_trans_ch
import copy
from Bio import SeqIO
import argparse
from PIL import Image
from my_error import *


def draw_Repetitive_sequence(projectid, outdir, inputfile):
    path = os.path.dirname(os.path.realpath(__file__))
    cmd = f'perl /tmp/repeat_disc/repeat_disc.pl {projectid} {outdir} {inputfile} " 1-10 2-6 3-5 4-5 5-5 6-5 " " 2 7 7 80 10 50 500 -f -d -m " " -f -p -h 3 -l 30 "'
    return os.system(cmd)


def rm_dir1(outdir):
    cmd = f"rm -r {os.path.join(outdir, 'trans_cis_file')}"
    return os.system(cmd)


def rm_dir2(outdir):
    cmd = f"rm -r {os.path.join(outdir, 'trans_cis_picture')}"
    return os.system(cmd)


def main(inputfile, outdir, projectid, type, interval_ratio, dis_R_Threshold, ratio_cis_trans, draw_exon, draw_shadow=True
         , draw_three_repeat=True, draw_repeat_region=None, draw_misc_feature=None, draw_ncRNA=None,
         draw_3UTR=None, draw_5UTR=None,
         draw_NADH=True, draw_SDH=True, draw_COB=True, draw_COX=True, draw_ATP=True, draw_CCB=True, draw_RPO=True,
         draw_RPS=True, draw_RPL=True, draw_CLP=True, draw_ORF=True, draw_TRN=True, draw_RRN=True, draw_others=True,
         draw_ch_PAS=True, draw_ch_PSB=True, draw_ch_PET=True, draw_ch_ATP=True, draw_ch_NDH=True,
         draw_ch_RBC=True, draw_ch_RPO=True, draw_ch_RSP=True, draw_ch_RPL=True, draw_ch_CLP=True,
         draw_ch_YCF=True, draw_ch_TRN=True, draw_ch_RRN=True, draw_ch_others=True, out_type='png'
         ):
    print(1)
    # with open('/home/xyzhang/pmgdraw/PMGview/media/downloads/1.txt', 'w') as f:
    #     f.write('123')
    trans_cis_picture_type = []
    # dna1
    for i in range(len(inputfile)):
        record = SeqIO.read(inputfile[i], "genbank")  # 读取gb文件
        if 'accessions' not in record.annotations.keys():
            title = f'contig_{i + 1}'
        else:
            title = record.annotations['accessions'][0]

        # 测试一下没有标注的情况
        if 'topology' not in record.annotations.keys():
            trans_cis_picture_type.append('circular')
        else:
            trans_cis_picture_type.append(record.annotations['topology'])
        if not os.path.isdir(os.path.join(outdir, title)):
            os.makedirs(os.path.join(outdir, title))
        file = open(os.path.join(outdir, title, title + '.fasta'), 'w')
        file.write('>' + title + '\n' + str(record.seq))
        file.close()
        record.annotations['organism'] = record.annotations['organism'].replace(' ', '\n            ', 1)
        # record2 = dna1.shorten_distance2(copy.deepcopy(record), ratio)
        # with open(os.path.join(outdir, title, f'{title}_{ratio}.gb'), 'w') as f:
        #     dna1.myGenBankWriter(f).write_file([record2])

        draw_Repetitive_sequence(projectid, os.path.join(outdir, title, 'Repetitive_sequence'),
                                 os.path.join(outdir, title, title + '.fasta'))

        # 散在重复序列
        dis_R = open(os.path.join(outdir, title, 'Repetitive_sequence', f'{projectid}_vmatch.txt'), 'r')  # 读取散在重复序列文件
        dis_file = dis_R.readlines()
        dis_file_new = open(os.path.join(outdir, title, title + '_dis_R.txt'), 'w')
        for j in range(1, len(dis_file)):
            if dis_file[j].split()[1] != dis_file[j].split()[4]:
                if dis_file[j].split()[2] == 'D':
                    dis_file_new.write(
                        dis_file[j].split()[0] + '\t' + dis_file[j].split()[1] + '\t' + '1' + '\t' +
                        dis_file[j].split()[4] + '\n')
                if dis_file[j].split()[2] == 'P':
                    dis_file_new.write(
                        dis_file[j].split()[0] + '\t' + dis_file[j].split()[1] + '\t' + '-1' + '\t' +
                        dis_file[j].split()[4] + '\n')
        dis_file_new.close()

        # 微卫星重复序列
        short_TR = open(os.path.join(outdir, title, 'Repetitive_sequence', f'{projectid}.fas.misa'), 'r')  # 读取微卫星重复序列文件
        short_file = short_TR.readlines()
        short_file_new = open(os.path.join(outdir, title, title + '_short_TR.txt'), 'w')
        for j in range(1, len(short_file)):
            if short_file[j].split()[5] != short_file[j].split()[6]:
                short_file_new.write(short_file[j].split()[5] + '\t' + short_file[j].split()[6] + '\n')
        short_file_new.close()

        # 串联重复序列
        long_TR = open(os.path.join(outdir, title, 'Repetitive_sequence', f'{projectid}.fas.2.7.7.80.10.50'
                                                                          f'.500.dat'),
                       'r')  # 读取串联重复序列文件
        long_file = long_TR.readlines()
        long_file_new = open(os.path.join(outdir, title, title + '_long_TR.txt'), 'w')
        for j in range(15, len(long_file)):
            if long_file[j].split()[0] != long_file[j].split()[1]:
                long_file_new.write(long_file[j].split()[0] + '\t' + long_file[j].split()[1] + '\n')
        long_file_new.close()

        if type == 'Mitochondrial':
            if trans_cis_picture_type[i] == 'circular':
                # draw circle trans_cis_picture
                circle_draw_new.my_circle_draw(inputfile=inputfile[i],
                                               short_TR=os.path.join(outdir, title, title + '_short_TR.txt'),
                                               long_TR=os.path.join(outdir, title, title + '_long_TR.txt'),
                                               dis_R=os.path.join(outdir, title, title + '_dis_R.txt'),
                                               outputfile=os.path.join(outdir, title + '_circular.' + out_type),
                                               dis_R_Threshold=dis_R_Threshold,
                                               interval_ratio=interval_ratio,
                                               draw_exon=draw_exon,
                                               draw_shadow=draw_shadow,
                                               draw_NADH=draw_NADH, draw_SDH=draw_SDH, draw_COB=draw_COB, draw_COX=draw_COX,
                                               draw_ATP=draw_ATP,
                                               draw_CCB=draw_CCB, draw_RPO=draw_RPO, draw_RPS=draw_RPS, draw_RPL=draw_RPL,
                                               draw_CLP=draw_CLP, draw_ORF=draw_ORF,
                                               draw_TRN=draw_TRN, draw_RRN=draw_RRN, draw_repeat_region=draw_repeat_region,
                                               draw_misc_feature=draw_misc_feature, draw_ncRNA=draw_ncRNA,
                                               draw_3UTR=draw_3UTR, draw_5UTR=draw_5UTR, draw_others=draw_others,
                                               draw_three_repeat=draw_three_repeat
                                               )
            elif trans_cis_picture_type[i] == 'linear':
                # draw linear trans_cis_picture
                linear_draw_new.my_linear_draw(inputfile=inputfile[i],
                                               short_TR=os.path.join(outdir, title, title + '_short_TR.txt'),
                                               long_TR=os.path.join(outdir, title, title + '_long_TR.txt'),
                                               dis_R=os.path.join(outdir, title, title + '_dis_R.txt'),
                                               outputfile=os.path.join(outdir, title + '_linear.' + out_type),
                                               dis_R_Threshold=dis_R_Threshold,
                                               interval_ratio=interval_ratio,
                                               draw_exon=draw_exon,
                                               draw_shadow=draw_shadow,
                                               draw_NADH=draw_NADH, draw_SDH=draw_SDH, draw_COB=draw_COB, draw_COX=draw_COX,
                                               draw_ATP=draw_ATP,
                                               draw_CCB=draw_CCB, draw_RPO=draw_RPO, draw_RPS=draw_RPS, draw_RPL=draw_RPL,
                                               draw_CLP=draw_CLP, draw_ORF=draw_ORF,
                                               draw_TRN=draw_TRN, draw_RRN=draw_RRN, draw_repeat_region=draw_repeat_region,
                                               draw_misc_feature=draw_misc_feature, draw_ncRNA=draw_ncRNA,
                                               draw_3UTR=draw_3UTR, draw_5UTR=draw_5UTR, draw_others=draw_others,
                                               draw_three_repeat=draw_three_repeat)

        elif type == 'Plastid':
            if trans_cis_picture_type[i] == 'circular':
                # draw circle trans_cis_picture
                circle_draw_ch_new.my_circle_draw(inputfile=inputfile[i],
                                                  short_TR=os.path.join(outdir, title, title + '_short_TR.txt'),
                                                  long_TR=os.path.join(outdir, title, title + '_long_TR.txt'),
                                                  dis_R=os.path.join(outdir, title, title + '_dis_R.txt'),
                                                  outputfile=os.path.join(outdir, title + '_circular.' + out_type),
                                                  dis_R_Threshold=dis_R_Threshold,
                                                  interval_ratio=interval_ratio,
                                                  draw_exon=draw_exon,
                                                  draw_shadow=draw_shadow, draw_repeat_region=draw_repeat_region, draw_misc_feature=draw_misc_feature,
                                                  draw_ncRNA=draw_ncRNA,
                                                  draw_3UTR=draw_3UTR, draw_5UTR=draw_5UTR, draw_three_repeat=draw_three_repeat,
                                                  draw_ch_PAS=draw_ch_PAS,
                                                  draw_ch_PSB=draw_ch_PSB, draw_ch_PET=draw_ch_PET, draw_ch_ATP=draw_ch_ATP,
                                                  draw_ch_NDH=draw_ch_NDH, draw_ch_RBC=draw_ch_RBC,
                                                  draw_ch_RPO=draw_ch_RPO, draw_ch_RSP=draw_ch_RSP, draw_ch_RPL=draw_ch_RPL,
                                                  draw_ch_CLP=draw_ch_CLP, draw_ch_YCF=draw_ch_YCF,
                                                  draw_ch_TRN=draw_ch_TRN, draw_ch_RRN=draw_ch_RRN, draw_ch_others=draw_ch_others)
            elif trans_cis_picture_type[i] == 'linear':
                # draw linear trans_cis_picture
                linear_draw_ch_new.my_linear_draw(inputfile=inputfile[i],
                                                  short_TR=os.path.join(outdir, title, title + '_short_TR.txt'),
                                                  long_TR=os.path.join(outdir, title, title + '_long_TR.txt'),
                                                  dis_R=os.path.join(outdir, title, title + '_dis_R.txt'),
                                                  outputfile=os.path.join(outdir, title + '_linear.' + out_type),
                                                  dis_R_Threshold=dis_R_Threshold,
                                                  interval_ratio=interval_ratio,
                                                  draw_exon=draw_exon,
                                                  draw_shadow=draw_shadow, draw_repeat_region=draw_repeat_region, draw_misc_feature=draw_misc_feature,
                                                  draw_ncRNA=draw_ncRNA,
                                                  draw_3UTR=draw_3UTR, draw_5UTR=draw_5UTR, draw_three_repeat=draw_three_repeat,
                                                  draw_ch_PAS=draw_ch_PAS,
                                                  draw_ch_PSB=draw_ch_PSB, draw_ch_PET=draw_ch_PET, draw_ch_ATP=draw_ch_ATP,
                                                  draw_ch_NDH=draw_ch_NDH, draw_ch_RBC=draw_ch_RBC,
                                                  draw_ch_RPO=draw_ch_RPO, draw_ch_RSP=draw_ch_RSP, draw_ch_RPL=draw_ch_RPL,
                                                  draw_ch_CLP=draw_ch_CLP, draw_ch_YCF=draw_ch_YCF,
                                                  draw_ch_TRN=draw_ch_TRN, draw_ch_RRN=draw_ch_RRN, draw_ch_others=draw_ch_others)
    dna3_ratio = 0

    # dna2
    print(2)
    trans_cis_out = []
    for i in range(len(inputfile)):
        if not os.path.isdir(os.path.join(outdir, 'trans_cis_file', 'trans_cis', str(i + 1), 'trans')):
            os.makedirs(os.path.join(outdir, 'trans_cis_file', 'trans_cis', str(i + 1), 'trans'))
        if not os.path.isdir(os.path.join(outdir, 'trans_cis_file', 'trans_cis', str(i + 1), 'cis')):
            os.makedirs(os.path.join(outdir, 'trans_cis_file', 'trans_cis', str(i + 1), 'cis'))
        trans_cis_out.append(os.path.join(outdir, 'trans_cis_file', 'trans_cis', str(i + 1)))
    records = [SeqIO.read(file_name, "genbank") for file_name in inputfile]
    for i in range(len(records)):
        records[i].annotations['organism'] = records[i].annotations['organism'].replace(' ', '\n            ', 1)
    dna2.shorten_distance(copy.deepcopy(records), ratio_cis_trans, trans_cis_out)
    for i in range(len(inputfile)):
        if not os.path.isdir(os.path.join(outdir, 'trans_cis_picture', 'trans_cis', str(i + 1), 'trans')):
            os.makedirs(os.path.join(outdir, 'trans_cis_picture', 'trans_cis', str(i + 1), 'trans'))
        trans_cis_file_names = os.listdir(os.path.join(outdir, 'trans_cis_file', 'trans_cis', str(i + 1), 'trans'))
        for trans_cis_file_name in trans_cis_file_names:
            if trans_cis_file_name[-3:] == '.gb':
                new_drawgenemap.mydrawgenemap(
                    os.path.join(outdir, 'trans_cis_file', 'trans_cis', str(i + 1), 'trans', trans_cis_file_name),
                    os.path.join(outdir, 'trans_cis_picture', 'trans_cis', str(i + 1), 'trans',
                                 trans_cis_file_name[:-3] + '.jpg'))

        if not os.path.isdir(os.path.join(outdir, 'trans_cis_picture', 'trans_cis', str(i + 1), 'cis')):
            os.makedirs(os.path.join(outdir, 'trans_cis_picture', 'trans_cis', str(i + 1), 'cis'))
        trans_cis_file_names = os.listdir(os.path.join(outdir, 'trans_cis_file', 'trans_cis', str(i + 1), 'cis'))
        for trans_cis_file_name in trans_cis_file_names:
            if trans_cis_file_name[-3:] == '.gb':
                new_drawgenemap.mydrawgenemap(
                    os.path.join(outdir, 'trans_cis_file', 'trans_cis', str(i + 1), 'cis', trans_cis_file_name),
                    os.path.join(outdir, 'trans_cis_picture', 'trans_cis', str(i + 1), 'cis',
                                 trans_cis_file_name[:-3] + '.jpg'))

    # mRNA
    if not os.path.isdir(os.path.join(outdir, 'trans_cis_file', 'mRNA', 'trans')):
        os.makedirs(os.path.join(outdir, 'trans_cis_file', 'mRNA', 'trans'))
    if not os.path.isdir(os.path.join(outdir, 'trans_cis_file', 'mRNA', 'cis')):
        os.makedirs(os.path.join(outdir, 'trans_cis_file', 'mRNA', 'cis'))
    records = [SeqIO.read(file_name, "genbank") for file_name in inputfile]
    for i in range(len(records)):
        records[i].annotations['organism'] = records[i].annotations['organism'].replace(' ', '\n            ', 1)
    dna3.shorten_distance(records, dna3_ratio, os.path.join(outdir, 'trans_cis_file', 'mRNA'))

    if not os.path.isdir(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'trans')):
        os.makedirs(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'trans'))
    file_names = os.listdir(os.path.join(outdir, 'trans_cis_file', 'mRNA', 'trans'))
    for file_name in file_names:
        if file_name[-3:] == '.gb':
            new_drawgenemap.mydrawgenemap(os.path.join(outdir, 'trans_cis_file', 'mRNA', 'trans', file_name),
                                          os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'trans',
                                                       file_name[:-3] + '.jpg'))
    if not os.path.isdir(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'cis')):
        os.makedirs(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'cis'))
    trans_cis_file_names = os.listdir(os.path.join(outdir, 'trans_cis_file', 'mRNA', 'cis'))
    for file_name in trans_cis_file_names:
        if file_name[-3:] == '.gb':
            print(file_name)
            if file_name == 'nad2  .gb':
                print(1)
            new_drawgenemap.mydrawgenemap(os.path.join(outdir, 'trans_cis_file', 'mRNA', 'cis', file_name),
                                          os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'cis',
                                                       file_name[:-3] + '.jpg'))

    record = SeqIO.read(inputfile[0], "genbank")
    if 'organism' not in record.annotations.keys():
        title = 'species'
    else:
        title = record.annotations['organism']

    CHRs = []
    for i in range(len(records)):
        if 'accessions' not in record.annotations.keys():
            CHR = f'contig_{i + 1}'
        else:
            CHR = records[i].annotations['accessions'][0]
        CHRs.append(CHR)
    if type == 'Mitochondrial':
        # trans
        if len(os.listdir(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'trans'))) > 0:
            header_img = my_cis_and_trans.obtain_gene_names(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'trans'),
                                                            'trans',
                                                            f'{title}')
            my_cis_and_trans.load_input2_fig(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'trans'))
            for i in range(len(inputfile)):
                my_cis_and_trans.load_input1_fig(
                    os.path.join(outdir, 'trans_cis_picture', 'trans_cis', str(i + 1), 'trans'),
                    CHR=CHRs[i])
            img = Image.fromarray(my_cis_and_trans.splice_genes(header_img))
            img.save(os.path.join(outdir, 'trans_splicing_genes.' + out_type))

        # cis
        if len(os.listdir(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'cis'))) > 0:
            header_img = my_cis_and_trans.obtain_gene_names(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'cis'),
                                                            'cis',
                                                            f'{title}')
            my_cis_and_trans.load_input2_fig(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'cis'))
            for i in range(len(inputfile)):
                my_cis_and_trans.load_input1_fig(
                    os.path.join(outdir, 'trans_cis_picture', 'trans_cis', str(i + 1), 'cis'),
                    CHR=CHRs[i])
            img = Image.fromarray(my_cis_and_trans.splice_genes(header_img))

            img.save(os.path.join(outdir, 'cis_splicing_genes.' + out_type))
    elif type == 'Plastid':
        # trans
        print(os.listdir(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'trans')))
        if len(os.listdir(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'trans'))) > 0:
            header_img = my_cis_and_trans_ch.obtain_gene_names(
                os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'trans'),
                'trans',
                f'{title}')
            my_cis_and_trans_ch.load_input2_fig(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'trans'))
            for i in range(len(inputfile)):
                my_cis_and_trans_ch.load_input1_fig(
                    os.path.join(outdir, 'trans_cis_picture', 'trans_cis', str(i + 1), 'trans'),
                    CHR=CHRs[i])
            img = Image.fromarray(my_cis_and_trans_ch.splice_genes(header_img))
            img.save(os.path.join(outdir, 'trans_splicing_genes.' + out_type))

        # cis
        if len(os.listdir(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'cis'))) > 0:
            header_img = my_cis_and_trans_ch.obtain_gene_names(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'cis'),
                                                               'cis',
                                                               f'{title}')
            my_cis_and_trans_ch.load_input2_fig(os.path.join(outdir, 'trans_cis_picture', 'mRNA', 'cis'))
            for i in range(len(inputfile)):
                my_cis_and_trans_ch.load_input1_fig(
                    os.path.join(outdir, 'trans_cis_picture', 'trans_cis', str(i + 1), 'cis'),
                    CHR=CHRs[i])
            img = Image.fromarray(my_cis_and_trans_ch.splice_genes(header_img))

            img.save(os.path.join(outdir, 'cis_splicing_genes.' + out_type))

    rm_dir1(outdir)
    rm_dir2(outdir)


if __name__ == "__main__":
    dna3_ratio = 0


    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise TypeError('Boolean value expected.')


    parser = argparse.ArgumentParser()
    parser.add_argument('--inputfile', help='the genbank file to input', type=str, default='/mnt/biodata02/home/xyzhang/pmgdraw-test/mt_transform/NC_026975.1.gb')
    parser.add_argument('--outdir', help='the path to save result', default='./123')
    parser.add_argument('--projectid', help='the path to save result', type=int, default=123)
    parser.add_argument('--type', help='Mitochondrial or Plastid', type=str, default='Mitochondrial')
    parser.add_argument('--out_type', help='Mitochondrial or Plastid', type=str, default='png')
    parser.add_argument('--interval_ratio', help='the ratio of interval to be kept,0-1', type=float, default=1)
    parser.add_argument('--dis-R-Threshold',
                        help='Scattered repeats greater than this threshold will be drawn in the outermost circle',
                        type=int, default=5000)
    parser.add_argument('--ratio-cis-trans', help='the ratio of interval to be kept,0-1', type=float, default=0.02)
    parser.add_argument('--draw_exon', help='the ratio of interval to be kept,0-1', type=str2bool, default=True)
    parser.add_argument('--draw_shadow', help='the ratio of interval to be kept,0-1', type=str2bool, default=True)
    args = parser.parse_args()
    args.inputfile = args.inputfile.split(',')
    args.trans_cis_picture_type = []

    try:
        main(args.inputfile, args.outdir, args.projectid, args.type, args.interval_ratio, args.dis_R_Threshold,
             args.ratio_cis_trans, args.draw_exon, args.draw_shadow)
    except noGeneNameError as e:
        print(1)
