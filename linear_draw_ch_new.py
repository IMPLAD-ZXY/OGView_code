import numpy as np
import Bio
from PIL import Image, ImageDraw, ImageFont
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.pyplot import Line2D
from my_error import *


colors = [(42 / 255, 99 / 255, 50 / 255), (76 / 255, 136 / 255, 5 / 255), (127 / 255, 153 / 255, 44 / 255),
          (159 / 255, 187 / 255, 61 / 255),
          (254 / 255, 238 / 255, 80 / 255), (77 / 255, 158 / 255, 63 / 255), (174 / 255, 45 / 255, 41 / 255),
          (214 / 255, 173 / 255, 124 / 255),
          (156 / 255, 122 / 255, 75 / 255), (217 / 255, 102 / 255, 45 / 255), (113 / 255, 184 / 255, 169 / 255),
          (23 / 255, 44 / 255, 127 / 255),
          (209 / 255, 56 / 255, 42 / 255), (125 / 255, 125 / 255, 125 / 255)]

def colors_index(name):
    if 'psa' in name or 'PSA' in name:
        color = 0

    elif 'psb' in name or 'PSB' in name:
        color = 1

    elif 'pet' in name or 'PET' in name:
        color = 2

    elif 'atp' in name or 'ATP' in name:
        color = 3

    elif 'ndh' in name or 'NDH' in name:
        color = 4

    elif 'rbc' in name or 'RBC' in name:
        color = 5
    elif 'rpo' in name or 'RPO' in name:
        color = 6
    elif 'rps' in name or 'RPS' in name:
        color = 7
    elif 'rpl' in name or 'RPL' in name:
        color = 8
    elif 'clp' in name or 'mat' in name or 'CLP' in name or 'MAT' in name:
        color = 9
    elif 'ycf' in name or 'YCF' in name:
        color = 10
    elif 'trn' in name or 'TRN' in name:
        color = 11
    elif 'rrn' in name or 'RRN' in name:
        color = 12
    else:
        color = 13
    return color

class genes_dict:
    def __init__(self):
        self.gene_names = []
        self.genes_info = {}
        self.gene_locations = []
        self.gene_location_bias = []
        self.length = None

    def add_gene(self, gene_name):
        # 如果gene_name不是self.genes_info的key，则添加
        if gene_name not in self.gene_names:
            self.gene_names.append(gene_name)
            self.genes_info[gene_name] = {'gene': [], 'PCGs': [], 'exon': []}

    def add_feature(self, feature):
        # 将gene、PCGs、exon类型的feature添加到genes_info中对应gene_name的字典中
        if feature.type == 'source':
            self.length = int(feature.location.end)
        if feature.type not in ['CDS', 'tRNA', 'rRNA', 'exon', 'gene']:
            return
        # assert len(feature.qualifiers['gene']) > 0, f'{feature} has no gene!!!'
        if 'gene' not in feature.qualifiers.keys() or len(feature.qualifiers['gene']) == 0:
            raise noGeneNameError(feature)
        gene_name = feature.qualifiers['gene'][0]
        self.add_gene(gene_name)
        if feature.type in ['CDS', 'tRNA', 'rRNA']:
            self.genes_info[gene_name]['PCGs'].append(feature)
        elif feature.type == 'exon':
            self.genes_info[gene_name]['exon'].append(feature)
        elif feature.type == 'gene':
            self.genes_info[gene_name]['gene'].append(feature)

    def compute_location_bias(self, interval_ratio):
        """
        :param interval_ratio:
        :return:
        """
        # 对所有gene类型feature的位置，按照起始位置从小到大排序，为后面调整间距做准备
        locations = []
        for gene_name in self.gene_names:
            for feature in self.genes_info[gene_name]['gene']:
                assert int(feature.location.end) <= self.length
                locations.append((int(feature.location.start), int(feature.location.end), gene_name, feature.strand))
        locations.sort(key=lambda location: location[0]-location[1]/(self.length+1))

        if len(locations) == 0:
            raise noAnnotatedGenomeError()

        self.gene_locations = []

        genes_in_location = []
        for i in range(len(locations)-1):
            gene_name = locations[i][2]
            if locations[i][1] > locations[i + 1][0] and len(self.genes_info[gene_name]['exon']) > 1:  # 如果跨过下一个gene
                if gene_name not in genes_in_location:
                    for j, feature in enumerate(self.genes_info[gene_name]['exon']):
                        self.gene_locations.append((int(feature.location.start), int(feature.location.end), gene_name+f'(exon{feature.qualifiers["number"][0]})', feature.strand))
                    genes_in_location.append(gene_name)
            else:
                self.gene_locations.append(locations[i])
        self.gene_locations.append(locations[-1])  # 刚刚没有加入最后一个gene的位置
        self.gene_locations.sort(key=lambda location: location[0]-location[1]/(self.length+1))  #在起点相同的情况下，让终点靠后的基因排在前面

        # 计算向前平移的距离
        last_end = -1
        bias = 0
        for i in range(len(self.gene_locations)):
            if self.gene_locations[i][0] >= last_end:
                bias += int((self.gene_locations[i][0] - last_end) * (1-interval_ratio))
            # else:
            #     assert False
            last_end = max(self.gene_locations[i][1], last_end)
            self.gene_location_bias.append(bias)

def my_linear_draw(inputfile, short_TR, long_TR, dis_R, outputfile, dis_R_Threshold, interval_ratio, draw_exon=True,
                   draw_shadow=True, draw_repeat_region=None, draw_misc_feature=None, draw_ncRNA=None,
                   draw_3UTR=None, draw_5UTR=None, draw_three_repeat=True, draw_ch_PAS=True,
                   draw_ch_PSB=True, draw_ch_PET=True, draw_ch_ATP=True, draw_ch_NDH=True, draw_ch_RBC=True,
                   draw_ch_RPO=True, draw_ch_RSP=True, draw_ch_RPL=True, draw_ch_CLP=True, draw_ch_YCF=True,
                   draw_ch_TRN=True, draw_ch_RRN=True, draw_ch_others=True):
    feature_types_to_shorten = ['CDS', 'tRNA', 'rRNA']

    # 读取 inputfile
    record = SeqIO.read(inputfile, "genbank")
    genes = genes_dict()

    for feature in record.features:
        genes.add_feature(feature)

    genes.compute_location_bias(interval_ratio)

    old_len = len(record.seq)
    new_len = int(old_len - genes.gene_locations[-1][1]) * interval_ratio + genes.gene_locations[-1][1]

    def remove_overlap(data):
        """
        :param data: data[i] (第i个文字的调整后的中心，第i个文字的原中心，基因名，正反链)
        :return:
        """
        index2cluster_old = np.arange(len(data), dtype=int)
        index2cluster_new = index2cluster_old.copy()
        while True:

            for i in range(len(data) - 1):
                if data[i + 1][0] - data[i][0] < text_interval:
                    new_cluster = index2cluster_new[i]
                    old_cluster = index2cluster_old[i + 1]
                    index2cluster_new[index2cluster_old == old_cluster] = new_cluster
                    data[i + 1][0] = data[i][0] + text_interval

            if (index2cluster_new == index2cluster_old).all():
                break

            clusters = []
            for i in range(len(data)):
                cluster = index2cluster_new[i]
                if cluster not in clusters:
                    clusters.append(cluster)
                    center = (data[int(np.where(index2cluster_new == cluster)[0][0])][0] +
                              data[int(np.where(index2cluster_new == cluster)[0][-1])][0]) / 2

                    cluster_width = len(np.where(index2cluster_new == cluster)[0]) * text_interval
                    target_center = min(long - cluster_width / 2, max(cluster_width / 2, (
                            data[int(np.where(index2cluster_new == cluster)[0][0])][1] +
                            data[int(np.where(index2cluster_new == cluster)[0][-1])][1]) / 2))
                    for j in np.where(index2cluster_new == cluster)[0]:
                        data[j][0] += target_center - center
            index2cluster_old = index2cluster_new.copy()
        return data

    if 'organism' not in record.annotations.keys():
        title = f'species_name'
    else:
        title = record.annotations['organism']

    if draw_three_repeat and interval_ratio == 1:
        line_height = 0.7
    else:
        line_height = 0
    long = 5
    height = 5
    blank = 0
    line_width = 1.5
    bar_width = 0.1
    bar_height = 0.15
    label_size = 4
    text_line_width = 0.4
    text_interval = 0.05
    text_height = 0.2
    scale_height = 0.12
    title_size = 10
    mitochondrial_size = 8
    line_long = long - 2 * blank
    ratio = line_long / new_len
    ratio_old = line_long / old_len
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(111)
    ax.set_xlim(0, long)  # 宽
    ax.set_ylim(0, height)  # 高

    # 上方中轴线
    line1 = Line2D((blank, long - blank), (height / 2 + line_height, height / 2 + line_height), color='black',
                   linewidth=line_width)  # 起点x轴，终点x轴，起点y轴，终点y轴
    ax.add_line(line1)
    # 下方中轴线
    line2 = Line2D((blank, long - blank), (height / 2 - line_height, height / 2 - line_height), color='black',
                   linewidth=line_width)  # 起点x轴，终点x轴，起点y轴，终点y轴
    ax.add_line(line2)

    if draw_three_repeat and interval_ratio == 1:
        # 刻度实线
        line2 = Line2D((blank, long - blank),
                       (height / 2 - line_height + scale_height, height / 2 - line_height + scale_height),
                       color='black', linewidth=0.7)  # 起点x轴，终点x轴，起点y轴，终点y轴
        ax.add_line(line2)
        # GC含量实线
        line2 = Line2D((blank, long - blank),
                       (height / 2 - line_height + bar_height / 2 + scale_height,
                        height / 2 - line_height + bar_height / 2 + scale_height),
                       color='black', linewidth=0.4)  # 起点x轴，终点x轴，起点y轴，终点y轴
        ax.add_line(line2)
        # GC含量虚线
        line2 = Line2D((blank, long - blank),
                       (height / 2 - line_height + bar_height + scale_height,
                        height / 2 - line_height + bar_height + scale_height),
                       color='black', ls='--', linewidth=0.7)  # 起点x轴，终点x轴，起点y轴，终点y轴 需要加上线的粗度
        ax.add_line(line2)
        # 串联重复实线
        line2 = Line2D((blank, long - blank),
                       (height / 2 + line_height - 0.1,
                        height / 2 + line_height - 0.1),
                       color='black', linewidth=0.7)  # 起点x轴，终点x轴，起点y轴，终点y轴
        ax.add_line(line2)
        # 微卫星重复实线
        line2 = Line2D((blank, long - blank),
                       (height / 2 + line_height - 0.2,
                        height / 2 + line_height - 0.2),
                       color='black', linewidth=0.7)  # 起点x轴，终点x轴，起点y轴，终点y轴
        ax.add_line(line2)
        # 微卫星重复实线
        line2 = Line2D((blank, long - blank),
                       (height / 2 + line_height - 0.3,
                        height / 2 + line_height - 0.3),
                       color='black', ls='--', linewidth=0.7)  # 起点x轴，终点x轴，起点y轴，终点y轴
        ax.add_line(line2)

    # 图注
    plt.gca().add_patch(plt.Rectangle(xy=(0, 0.62), width=0.05, height=0.05, fill=True, color=colors[0]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.002, 0.62), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.1, 0.64, 'photosystem I', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(plt.Rectangle(xy=(0, 0.54), width=0.05, height=0.05, fill=True, color=colors[1]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.002, 0.54), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.1, 0.56, 'photosystem II', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(plt.Rectangle(xy=(0, 0.46), width=0.05, height=0.05, fill=True, color=colors[2]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.002, 0.46), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.1, 0.48, 'cytochrome b/f complex', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(plt.Rectangle(xy=(0, 0.38), width=0.05, height=0.05, fill=True, color=colors[3]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.002, 0.38), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.1, 0.40, 'ATP synthesis', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(plt.Rectangle(xy=(0, 0.30), width=0.05, height=0.05, fill=True, color=colors[4]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.002, 0.30), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.1, 0.32, 'NADH dehydrogenase', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(plt.Rectangle(xy=(0, 0.22), width=0.05, height=0.05, fill=True, color=colors[5]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.002, 0.22), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.1, 0.24, 'Rubiscp large subunit', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(plt.Rectangle(xy=(0, 0.14), width=0.05, height=0.05, fill=True, color=colors[6]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.002, 0.14), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.1, 0.16, 'RNA polymerase', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.text(0, 0.08, 'The colored parabola in the center circle represents the Dispersed Repeats', ha='left',
             va='center', rotation=0, wrap=True,
             fontsize=label_size)

    plt.gca().add_patch(
        plt.Rectangle(xy=(1.4, 0.62), width=0.05, height=0.05, fill=True, color=colors[7]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(1.402, 0.62), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(1.5, 0.64, 'small ribosomal protein', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(
        plt.Rectangle(xy=(1.4, 0.54), width=0.05, height=0.05, fill=True, color=colors[8]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(1.402, 0.54), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(1.5, 0.56, 'large ribosomal protein', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(
        plt.Rectangle(xy=(1.4, 0.46), width=0.05, height=0.05, fill=True, color=colors[9]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(1.402, 0.46), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(1.5, 0.48, 'clpP,matK,infA', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(
        plt.Rectangle(xy=(1.4, 0.38), width=0.05, height=0.05, fill=True, color=colors[10]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(1.402, 0.38), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(1.5, 0.40, 'hypothetical reading frame', ha='left', va='center', rotation=0, wrap=True,
             fontsize=label_size)
    plt.gca().add_patch(
        plt.Rectangle(xy=(1.4, 0.30), width=0.05, height=0.05, fill=True, color=colors[11], alpha=0.6))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(1.402, 0.30), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(1.5, 0.32, 'transfer RNA', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(
        plt.Rectangle(xy=(1.4, 0.22), width=0.05, height=0.05, fill=True, color=colors[12], alpha=0.6))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(1.402, 0.22), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(1.5, 0.24, 'ribosomal RNA', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(
        plt.Rectangle(xy=(1.4, 0.14), width=0.05, height=0.05, fill=True, color=colors[13]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(1.402, 0.14), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(1.5, 0.16, 'other', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)

    plt.gca().add_patch(
        plt.Rectangle(xy=(2.8, 0.62), width=0.05, height=0.05, fill=True, color='Tomato'))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(2.802, 0.62), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(2.9, 0.64, 'Short Tandem Repeats', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(
        plt.Rectangle(xy=(2.8, 0.54), width=0.05, height=0.05, fill=True, color='BlueViolet', alpha=0.6))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(2.802, 0.54), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(2.9, 0.56, 'Long Tandem Repeats', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(
        plt.Rectangle(xy=(2.8, 0.46), width=0.05, height=0.05, fill=True, color='orange', alpha=0.6))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(2.802, 0.46), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(2.9, 0.48, 'Dispersed Repeats(direct matches)', ha='left', va='center', rotation=0, wrap=True,
             fontsize=label_size)
    plt.gca().add_patch(
        plt.Rectangle(xy=(2.8, 0.38), width=0.05, height=0.05, fill=True, color='green', alpha=0.6))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(2.802, 0.38), width=0.05, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(2.9, 0.40, 'Dispersed Repeats(palindromic matches)', ha='left', va='center', rotation=0, wrap=True,
             fontsize=label_size)
    plt.gca().add_patch(plt.Rectangle(xy=(2.8,0.30), width=0.05,height=0.05, fill=True, color=(184/255,211/255,143/255), alpha=0.6))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(2.802, 0.30), width=0.05, height=0.05, fill=False, edgecolor='black',linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(2.9, 0.32, 'Inverted Repeat(IR)', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(plt.Rectangle(xy=(2.8,0.22), width=0.05,height=0.05, fill=True, color=colors[10], alpha=0.4))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(2.802, 0.22), width=0.05, height=0.05, fill=False, edgecolor='black',linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(2.9, 0.24, 'Large Single Copy(LSR)', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)
    plt.gca().add_patch(plt.Rectangle(xy=(2.8,0.14), width=0.05,height=0.05, fill=True, color=(238/255, 220/255, 130/255), alpha=0.5))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(2.802, 0.14), width=0.05, height=0.05, fill=False, edgecolor='black',linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(2.9, 0.16, 'Small Single Copy(SSR)', ha='left', va='center', rotation=0, wrap=True, fontsize=label_size)


    plt.scatter(4.2, 0.64, s=10, edgecolors=(226 / 255, 0 / 255, 26 / 255), marker='o', alpha=0.8, linewidths=0.6, facecolors='none',)
    plt.text(4.4, 0.64, 'misc feature', ha='left', va='center', rotation=0, wrap=True,fontsize=label_size)
    plt.scatter(4.2, 0.56, s=10, edgecolors=(255 / 255, 105 / 255, 180 / 255), marker='*', alpha=0.8, linewidths=0.6,
                facecolors='none')
    plt.text(4.4, 0.56, 'ncRNA', ha='left', va='center', rotation=0, wrap=True,fontsize=label_size)
    plt.scatter(4.2, 0.48, s=10, edgecolors=(102 / 255, 205 / 255, 170 / 255), marker='v', alpha=0.8, linewidths=0.6,
                facecolors='none')
    plt.text(4.4, 0.48, 'repeat region', ha='left', va='center', rotation=0, wrap=True,fontsize=label_size)
    plt.scatter(4.2, 0.40, s=10, edgecolors=(151 / 255, 190 / 255, 13 / 255), marker='P', alpha=0.8, linewidths=0.6, facecolors='none')
    plt.text(4.4, 0.40, "5'UTR", ha='left', va='center', rotation=0, wrap=True,fontsize=label_size)
    plt.scatter(4.2, 0.32, s=10, edgecolors='brown', marker='d', alpha=0.8, linewidths=0.6, facecolors='none')
    plt.text(4.4, 0.32, "3'UTR", ha='left', va='center', rotation=0, wrap=True,fontsize=label_size)



    # 总的GC含量计算 表头标签
    G_C_con = record.seq.count('G')+record.seq.count('C')

    lens = str('{:,}'.format(old_len)) + ' bp  GC: ' + str("%.2f" % (100 * G_C_con / len(record.seq))) + '%'
    plt.text(long / 2, height - 0.2, title, ha='center', va='center', rotation=0,
             wrap=True, fontsize=title_size, fontweight='bold', style='italic')
    plt.text(long / 2, height - 0.4, 'Plastid Genome', ha='center', va='center', rotation=0,
             wrap=True, fontsize=mitochondrial_size)
    plt.text(long / 2, height - 0.6, lens, ha='center', va='center', rotation=0,
             wrap=True, fontsize=mitochondrial_size)

    if draw_three_repeat and interval_ratio == 1:
        # GC含量
        GC_content = []
        for i in range(int(old_len/100)):
            G_content = record.seq[100*i:(i+1)*100].count('G')
            C_content = record.seq[100*i:(i+1)*100].count('C')
            GC_content.append(G_content + C_content)

        # 画GC含量直方图
        width = (long - blank) / len(GC_content)
        for i in range(len(GC_content)):
            line1 = Line2D((i*width, i*width), (height / 2 - line_height + scale_height+bar_height/2, height / 2 - line_height + scale_height+bar_height/2+bar_height*(GC_content[i]-50)/100), color='orange',
                           linewidth=0.5)  # 起点x轴，终点x轴，起点y轴，终点y轴
            ax.add_line(line1)

        # 画刻度
        for i in range(12):
            scale_len = old_len / 12 * i * ratio_old
            scale_label = str(int(old_len / 12000 * i)) + 'kb'
            line2 = Line2D((scale_len, scale_len),
                           (height / 2 - line_height, height / 2 - line_height + 0.03), color='black',
                           linewidth=0.5)  # 起点x轴，终点x轴，起点y轴，终点y轴
            ax.add_line(line2)
            plt.text(scale_len, height / 2 - line_height + 0.06, scale_label, ha='center', va='center', rotation=0,
                     wrap=True, fontsize=4, color='blue')

    label_location1 = []
    label_location2 = []
    genes_flag = {}
    for gene_name in genes.gene_names:
        genes_flag[gene_name] = False
    for i, location in enumerate(genes.gene_locations):
        gene_name = location[2]
        if (draw_ch_PAS or colors_index(gene_name) != 0) \
                and (draw_ch_PSB or colors_index(gene_name) != 1) \
                and (draw_ch_PET or colors_index(gene_name) != 2) \
                and (draw_ch_ATP or colors_index(gene_name) != 3) \
                and (draw_ch_NDH or colors_index(gene_name) != 4) \
                and (draw_ch_RBC or colors_index(gene_name) != 5) \
                and (draw_ch_RPO or colors_index(gene_name) != 6) \
                and (draw_ch_RSP or colors_index(gene_name) != 7) \
                and (draw_ch_RPL or colors_index(gene_name) != 8) \
                and (draw_ch_CLP or colors_index(gene_name) != 9) \
                and (draw_ch_YCF or colors_index(gene_name) != 10) \
                and (draw_ch_TRN or colors_index(gene_name) != 11) \
                and (draw_ch_RRN or colors_index(gene_name) != 12) \
                and (draw_ch_others or colors_index(gene_name) != 13):
            if gene_name in genes.gene_names and len(genes.genes_info[gene_name]['exon']) > 0 and draw_exon:
                if not genes_flag[gene_name]:
                    for feature in genes.genes_info[gene_name]['exon']:
                        bias = 0
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        for j in range(len(genes.gene_locations)):
                            if genes.gene_locations[j][0] > start:
                                break
                            else:
                                bias = genes.gene_location_bias[j]

                        if feature.strand == 1:
                            plt.gca().add_patch(
                                plt.Rectangle(xy=(
                                ratio * (int(start) - bias) + blank, height / 2 + line_height),
                                              width=ratio * (int(end) - int(start)),
                                              height=bar_height, fill=True, color=colors[colors_index(gene_name)],
                                              linewidth=bar_width))  # xy是矩形的左下角坐标
                            plt.gca().add_patch(
                                plt.Rectangle(xy=(ratio * int(start) + blank, height / 2 + line_height),
                                              width=ratio * (int(end) - int(start)),
                                              height=bar_height, fill=False, linewidth=bar_width))  # xy是矩形的左下角坐标
                            center = ratio * (int(start) + int(end)) / 2 + blank  # 为什么要加0.02
                            label_location1.append([center, center, gene_name+f'(exon{int(feature.qualifiers["number"][0])})', 1])
                        else:
                            plt.gca().add_patch(
                                plt.Rectangle(xy=(
                                ratio * (int(start) - bias) + blank, height / 2 - line_height),
                                              width=ratio * (int(end) - int(start)),
                                              height=-bar_height, fill=True, color=colors[colors_index(gene_name)],
                                              linewidth=bar_width))  # xy是矩形的左下角坐标
                            plt.gca().add_patch(
                                plt.Rectangle(xy=(ratio * int(start) + blank, height / 2 - line_height),
                                              width=ratio * (int(end) - int(start)),
                                              height=-bar_height, fill=False, linewidth=bar_width))  # xy是矩形的左下角坐标
                            center = ratio * (int(start) + int(end)) / 2 + blank
                            label_location2.append([center + 0.02, center, gene_name+f'(exon{int(feature.qualifiers["number"][0])})', -1])
                            # plt.text(ratio * (location[0] + location[1]) / 2 + blank - 0.02, height / 2 - bar_height - 0.1, i,
                            #          ha='left', rotation=-90,
                            #          wrap=True, fontsize=label_size)
                    for feature in genes.genes_info[gene_name]['PCGs']:
                        if type(feature.location) != Bio.SeqFeature.CompoundLocation:
                            bias = 0
                            start = int(feature.location.start)
                            end = int(feature.location.end)
                            draw_flag = True
                            for exon_feature in genes.genes_info[gene_name]['exon']:
                                if start == int(exon_feature.location.start) and end == int(exon_feature.location.end):
                                    draw_flag = False
                                    break
                            if draw_flag:
                                for j in range(len(genes.gene_locations)):
                                    if genes.gene_locations[j][0] > start:
                                        break
                                    else:
                                        bias = genes.gene_location_bias[j]
                                if feature.strand == 1:
                                    plt.gca().add_patch(
                                        plt.Rectangle(xy=(
                                        ratio * (int(start) - bias) + blank, height / 2 + line_height),
                                                      width=ratio * (int(end) - int(start)),
                                                      height=bar_height, fill=True, color=colors[colors_index(gene_name)],
                                                      linewidth=bar_width))  # xy是矩形的左下角坐标
                                    plt.gca().add_patch(
                                        plt.Rectangle(xy=(ratio * int(start) + blank, height / 2 + line_height),
                                                      width=ratio * (int(end) - int(start)),
                                                      height=bar_height, fill=False, linewidth=bar_width))  # xy是矩形的左下角坐标
                                    center = ratio * (int(start) + int(end)) / 2 + blank  # 为什么要加0.02
                                    label_location1.append([center, center, gene_name, 1])
                                else:
                                    plt.gca().add_patch(
                                        plt.Rectangle(xy=(
                                        ratio * (int(start) - bias) + blank, height / 2 - line_height),
                                                      width=ratio * (int(end) - int(start)),
                                                      height=-bar_height, fill=True, color=colors[colors_index(gene_name)],
                                                      linewidth=bar_width))  # xy是矩形的左下角坐标
                                    plt.gca().add_patch(
                                        plt.Rectangle(xy=(ratio * int(start) + blank, height / 2 - line_height),
                                                      width=ratio * (int(end) - int(start)),
                                                      height=-bar_height, fill=False, linewidth=bar_width))  # xy是矩形的左下角坐标
                                    center = ratio * (int(start) + int(end)) / 2 + blank
                                    label_location2.append([center + 0.02, center, gene_name, -1])
                    genes_flag[gene_name] = True
            else:
                if location[3] == 1:
                    plt.gca().add_patch(
                        plt.Rectangle(xy=(ratio * (location[0] - genes.gene_location_bias[i]) + blank, height / 2 + line_height),
                                      width=ratio * (location[1] - location[0]),
                                      height=bar_height, fill=True, color=colors[colors_index(gene_name)],
                                      linewidth=bar_width))  # xy是矩形的左下角坐标
                    plt.gca().add_patch(
                        plt.Rectangle(xy=(ratio * location[0] + blank, height / 2 + line_height),
                                      width=ratio * (location[1] - location[0]),
                                      height=bar_height, fill=False, linewidth=bar_width))  # xy是矩形的左下角坐标
                    center = ratio * (location[0] + location[1]) / 2 + blank  # 为什么要加0.02
                    label_location1.append([center, center, gene_name, 1])
                else:
                    plt.gca().add_patch(
                        plt.Rectangle(xy=(ratio * (location[0] - genes.gene_location_bias[i]) + blank, height / 2 - line_height),
                                      width=ratio * (location[1] - location[0]),
                                      height=-bar_height, fill=True, color=colors[colors_index(gene_name)],
                                      linewidth=bar_width))  # xy是矩形的左下角坐标
                    plt.gca().add_patch(
                        plt.Rectangle(xy=(ratio * location[0] + blank, height / 2 - line_height),
                                      width=ratio * (location[1] - location[0]),
                                      height=-bar_height, fill=False, linewidth=bar_width))  # xy是矩形的左下角坐标
                    center = ratio * (location[0] + location[1]) / 2 + blank
                    label_location2.append([center + 0.02, center, gene_name, -1])
                    # plt.text(ratio * (location[0] + location[1]) / 2 + blank - 0.02, height / 2 - bar_height - 0.1, i,
                    #          ha='left', rotation=-90,
                    #          wrap=True, fontsize=label_size)


    # 排序
    for i in range(len(label_location1)):
        for index in range(len(label_location1) - 1):
            if label_location1[index][0] > label_location1[index + 1][0]:
                value = label_location1.pop(index)
                label_location1.insert(index + 1, value)

    for i in range(len(label_location2)):
        for index in range(len(label_location2) - 1):
            if label_location2[index][0] > label_location2[index + 1][0]:
                value = label_location2.pop(index)
                label_location2.insert(index + 1, value)

    label_location1 = remove_overlap(label_location1)
    for i in range(len(label_location1)):
        center, center_old, text, strand = label_location1[i]
        line = Line2D((center_old, center_old),
                      (height / 2 + line_height + bar_height, height / 2 + line_height + bar_height + text_height / 6),
                      color='black', linewidth=text_line_width)  # 起点x轴，终点x轴，起点y轴，终点y轴
        ax.add_line(line)
        line = Line2D((center_old, center),
                      (height / 2 + line_height + bar_height + text_height / 6,
                       height / 2 + line_height + bar_height + 5 * text_height / 6 - 0.02),
                      color='black', linewidth=text_line_width)  # 起点x轴，终点x轴，起点y轴，终点y轴
        ax.add_line(line)
        line = Line2D((center, center),
                      (height / 2 + line_height + bar_height + 5 * text_height / 6 - 0.02,
                       height / 2 + line_height + bar_height + text_height - 0.02),
                      color='black', linewidth=text_line_width)  # 起点x轴，终点x轴，起点y轴，终点y轴
        ax.add_line(line)
        plt.text(center + 0.025, height / 2 + line_height + bar_height + text_height, text, ha='left', rotation=90,
                 wrap=True, fontsize=label_size)  # 0.02是为了矫正写字时的坐标不是中心

    label_location2 = remove_overlap(label_location2)
    for i in range(len(label_location2)):
        center, center_old, text, strand = label_location2[i]
        line = Line2D((center_old, center_old), (
            height / 2 - line_height - bar_height - 0.003, height / 2 - line_height - bar_height - text_height / 6),
                      color='black', linewidth=text_line_width)  # 起点x轴，终点x轴，起点y轴，终点y轴
        ax.add_line(line)
        line = Line2D((center_old, center),
                      (height / 2 - line_height - bar_height - text_height / 6,
                       height / 2 - line_height - bar_height - 5 * text_height / 6 + 0.02),
                      color='black', linewidth=text_line_width)  # 起点x轴，终点x轴，起点y轴，终点y轴
        ax.add_line(line)
        line = Line2D((center, center),
                      (height / 2 - line_height - bar_height - 5 * text_height / 6 + 0.02,
                       height / 2 - line_height - bar_height - text_height + 0.02),
                      color='black', linewidth=text_line_width)  # 起点x轴，终点x轴，起点y轴，终点y轴
        ax.add_line(line)
        plt.text(center + 0.025, height / 2 - line_height - bar_height - text_height, text, ha='right', rotation=90,
                 wrap=True, fontsize=label_size)  # 0.02是为了矫正写字时的坐标不是中心

    if draw_three_repeat and interval_ratio == 1:
        # 画重复序列
        short_TR = open(short_TR, 'r')  # 读取微卫星重复序列文件
        short_location_old = short_TR.readlines()
        short_location = []
        for i in range(len(short_location_old)):
            short_location.append((int(short_location_old[i].replace('\n', '').split('\t')[0]),
                                   int(short_location_old[i].replace('\n', '').split('\t')[1])))

        long_TR = open(long_TR, 'r')  # 读取串联重复序列文件
        long_location_old = long_TR.readlines()
        long_location = []
        for i in range(len(long_location_old)):
            long_location.append((int(long_location_old[i].replace('\n', '').split('\t')[0]),
                                  int(long_location_old[i].replace('\n', '').split('\t')[1])))

        dis_R = open(dis_R, 'r') # 读取散在重复序列文件
        dis_location_old = dis_R.readlines()
        dis_location = []
        for i in range(len(dis_location_old)):
            dis_location.append((int(dis_location_old[i].replace('\n', '').split('\t')[0]),
                                 int(dis_location_old[i].replace('\n', '').split('\t')[1]),
                                 int(dis_location_old[i].replace('\n', '').split('\t')[2]),
                                 int(dis_location_old[i].replace('\n', '').split('\t')[3])))

        for i in short_location:
            plt.gca().add_patch(plt.Rectangle(xy=(ratio_old * i[0]+blank, height / 2 + line_height - 0.1), width=ratio_old * (i[1] - i[0]) / 2,
                                              height=0.1, fill=True, color='Tomato'))  # xy是矩形的左下角坐标

        for i in long_location:
            plt.gca().add_patch(plt.Rectangle(xy=(ratio_old * i[0]+blank, height / 2 + line_height - 0.2), width=ratio_old * (i[1] - i[0]) / 2,
                                              height=0.1, fill=True, color='BlueViolet'))  # xy是矩形的左下角坐标
            # plt.scatter(ratio_old * i[0]+blank, height / 2 + line_height, s=10, c='BlueViolet', marker='v', edgecolor='purple')


        for i in range(len(dis_location)):
            width1 = blank + ratio_old * (dis_location[i][1] + dis_location[i][0]/2)
            width2 = blank + ratio_old * (dis_location[i][3] + dis_location[i][0]/2)
            width4 = blank + ratio_old * dis_location[i][1]
            width5 = blank + ratio_old * (dis_location[i][1] + dis_location[i][0])
            width6 = blank + ratio_old * dis_location[i][3]
            width7 = blank + ratio_old * (dis_location[i][3] + dis_location[i][0])
            if dis_location[i][2] == 1:
                plt.gca().add_patch(
                    plt.Rectangle(xy=(width4, height / 2 + line_height - 0.3), width=ratio_old * dis_location[i][0],
                                height=0.1, fill=True, color='orange', alpha=0.6))  # xy是矩形的左下角坐标
                plt.gca().add_patch(
                    plt.Rectangle(xy=(width6, height / 2 + line_height - 0.3), width=ratio_old * dis_location[i][0],
                                height=0.1, fill=True, color='orange', alpha=0.6))  # xy是矩形的左下角坐标
            if dis_location[i][2] == -1:
                plt.gca().add_patch(
                    plt.Rectangle(xy=(width4, height / 2 + line_height - 0.3), width=ratio_old * dis_location[i][0],
                                height=0.1, fill=True, color='green', alpha=0.6))  # xy是矩形的左下角坐标
                plt.gca().add_patch(
                    plt.Rectangle(xy=(width6, height / 2 + line_height - 0.3), width=ratio_old * dis_location[i][0],
                                height=0.1, fill=True, color='green', alpha=0.6))  # xy是矩形的左下角坐标


            widths = []
            heights = []
            if width1 > width2:
                width1, width2 = width2, width1
            for x in np.arange(width1, width2, 0.01):
                height1 = height / 2 + line_height - 0.3
                height2 = height / 2 - line_height + bar_height + scale_height + 0.4 - 0.4 * (width2 - width1) / line_long
                a = 4 * (height1 - height2) / (width1 - width2) ** 2
                y = (a * x ** 2 - a * (width1 + width2) * x + height2 + a * (width1 + width2) ** 2 / 4)
                widths.append(x)
                heights.append(y)
            plt.plot(widths, heights, linewidth=0.5, alpha=0.8)

    if draw_shadow and draw_three_repeat and interval_ratio == 1:
        dis_length = []
        for i in range(len(dis_location)):
            dis_length.append(dis_location[i][0])
        max_index = np.argmax(np.array(dis_length))
        IR_A_start = dis_location[max_index][1]
        IR_A_end = dis_location[max_index][1] + dis_location[max_index][0]
        IR_B_start = dis_location[max_index][3]
        IR_B_end = dis_location[max_index][3] + dis_location[max_index][0]
        if IR_A_start > IR_B_start:
            IR_B_start, IR_A_start = IR_A_start, IR_B_start
        theta_IR_A = [ratio_old * IR_A_start, ratio_old * IR_A_end]
        theta_IR_B = [ratio_old * IR_B_start, ratio_old * IR_B_end]
        theta_LSC = [ratio_old * IR_B_end, ratio_old * IR_A_start + 2 * np.pi]
        theta_SSC = [ratio_old * IR_A_end, ratio_old * IR_B_start]
        if ratio_old * IR_A_start + line_long - ratio_old * IR_B_end < ratio_old * IR_B_start - ratio_old * IR_A_end:
            theta_SSC, theta_LSC = theta_LSC, theta_SSC

        h1 = height / 2 + line_height
        h2 = height / 2 - line_height
        a = [theta_IR_A[0], theta_IR_A[1], theta_IR_A[1], theta_IR_A[0]]
        b = [h1, h1, h2, h2]
        plt.fill(a, b, color=(184 / 255, 211 / 255, 143 / 255), alpha=0.5)  # 背景色

        a2 = [theta_IR_B[0], theta_IR_B[1], theta_IR_B[1], theta_IR_B[0]]
        plt.fill(a2, b, color=(184 / 255, 211 / 255, 143 / 255), alpha=0.5)  # 背景色

        if ratio_old * IR_A_start + line_long - ratio_old * IR_B_end < ratio_old * IR_B_start - ratio_old * IR_A_end:
            theta_SSC, theta_LSC = theta_LSC, theta_SSC
            a3 = [theta_SSC[0], line_long, line_long, theta_SSC[0]]
            plt.fill(a3, b, color=colors[10], alpha=0.3)  # 背景色
            a3_2 = [0, theta_SSC[1], theta_SSC[1], 0]
            plt.fill(a3_2, b, color=colors[10], alpha=0.3)  # 背景色
            a4 = [theta_LSC[0], theta_LSC[1], theta_LSC[1], theta_LSC[0]]
            plt.fill(a4, b, color=(238 / 255, 220 / 255, 130 / 255), alpha=0.4)  # 背景色
        else:
            a3 = [theta_LSC[0], line_long, line_long, theta_LSC[0]]
            plt.fill(a3, b, color=colors[10], alpha=0.3)  # 背景色
            a3_2 = [0, theta_LSC[1], theta_LSC[1], 0]
            plt.fill(a3_2, b, color=colors[10], alpha=0.3)  # 背景色
            a4 = [theta_SSC[0], theta_SSC[1], theta_SSC[1], theta_SSC[0]]
            plt.fill(a4, b, color=(238 / 255, 220 / 255, 130 / 255), alpha=0.4)  # 背景色

    for i in range(1, len(record.features)):
        if draw_repeat_region:
            if record.features[i].type == 'repeat_region':
                theta = ratio_old * (int(record.features[i].location.start) + int(record.features[i].location.end)) / 2
                plt.scatter(theta, height/2+line_height+0.07, s=10, edgecolors=(102 / 255, 205 / 255, 170 / 255), marker='v', alpha=0.8,
                            linewidths=0.6, facecolors='none')
        if draw_misc_feature:
            if record.features[i].type == 'misc_feature':
                theta = ratio_old * (int(record.features[i].location.start) + int(record.features[i].location.end)) / 2
                plt.scatter(theta, height/2+line_height, s=10, edgecolors=(226 / 255, 0 / 255, 26 / 255), marker='o', alpha=0.8, linewidths=0.6,
                            facecolors='none')
        if draw_ncRNA:
            if record.features[i].type == 'ncRNA':
                theta = ratio_old * (int(record.features[i].location.start) + int(record.features[i].location.end)) / 2
                plt.scatter(theta, height/2-line_height-0.07, s=10, edgecolors=(255 / 255, 105 / 255, 180 / 255), marker='*', alpha=0.8,
                            linewidths=0.6, facecolors='none')
        if draw_5UTR:
            if record.features[i].type == "5'UTR":
                theta = ratio_old * (int(record.features[i].location.start) + int(record.features[i].location.end)) / 2
                plt.scatter(theta, height/2-line_height-0.14, s=10, edgecolors=(151 / 255, 190 / 255, 13 / 255), marker='P', alpha=0.8, linewidths=0.6,
                            facecolors='none')
        if draw_3UTR:
            if record.features[i].type == "3'UTR":
                theta = ratio_old * (int(record.features[i].location.start) + int(record.features[i].location.end)) / 2
                plt.scatter(theta, height/2+line_height+0.14, s=13, edgecolors='brown', marker='d', alpha=0.8, linewidths=0.6, facecolors='none')

    plt.axis('off')
    plt.savefig(outputfile, bbox_inches='tight', pad_inches=0.5, dpi=200, transparent=True)
    # plt.show()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='the RNA file to edit', type=str,
                        default='./NC_037304.2.gb')
    parser.add_argument('--short_TR', help='the RNA file to edit', type=str,
                        default='./NC_037304_short_TR.txt')
    parser.add_argument('--long_TR', help='the RNA file to edit', type=str,
                        default='./NC_037304_long_TR.txt')
    parser.add_argument('--dispersed_R', help='the RNA file to edit', type=str,
                        default='./NC_037304_dis_R.txt')
    parser.add_argument('--output', help='the RNA file to edit', type=str,
                        default='./contig2_circle.png')
    parser.add_argument('--dis_R_Threshold',
                        help='Scattered repeats greater than this threshold will be drawn in the outermost circle',
                        type=int, default=5000)
    parser.add_argument('--interval_ratio', help='the ratio of interval to be kept,0-1', type=float, default=1)
    parser.add_argument('--draw_exon', help='the ratio of interval to be kept,0-1', type=bool, default=True)
    parser.add_argument('--draw_shadow', help='the ratio of interval to be kept,0-1', type=bool, default=True)
    parser.add_argument('--draw_three_repeat', help='the ratio of interval to be kept,0-1', type=bool, default=True)
    args = parser.parse_args()
    my_linear_draw(args.input, args.short_TR, args.long_TR, args.dispersed_R, args.output, args.dis_R_Threshold, args.interval_ratio, args.draw_exon, args.draw_shadow, draw_three_repeat=args.draw_three_repeat)
