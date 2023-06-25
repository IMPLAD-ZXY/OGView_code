import math
import time
import Bio
import matplotlib.patches
import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt
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


def draw_parabola(r1, r2, theta1, theta2):
    thetas3 = []
    r3_s = []
    if theta2 - theta1 == 0:
        print(1)
    a = (r2 * math.cos((theta2 - theta1) / 2) - r1) / (r2 * r2 * math.sin((theta2 - theta1) / 2) ** 2)
    theta = (theta2 + theta1) / 2
    for x in np.arange(-abs(r2 * math.sin((theta2 - theta1) / 2)), abs(r2 * math.sin((theta2 - theta1) / 2)), 0.05):
        y = a * x ** 2 + r1
        thetas3.append(np.arctan(x / y) + theta)
        r3_s.append(np.sqrt(x ** 2 + y ** 2))
    # thetas3.append(theta2)
    # r3_s.append(r2)
    return thetas3, r3_s


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
        locations.sort(key=lambda location: location[0] - location[1] / (self.length + 1))

        if len(locations) == 0:
            raise noAnnotatedGenomeError()

        self.gene_locations = []

        genes_in_location = []
        for i in range(len(locations) - 1):
            gene_name = locations[i][2]
            if locations[i][1] > locations[i + 1][0] and len(self.genes_info[gene_name]['exon']) > 1:  # 如果跨过下一个gene
                if gene_name not in genes_in_location:
                    for j, feature in enumerate(self.genes_info[gene_name]['exon']):
                        self.gene_locations.append((int(feature.location.start), int(feature.location.end),
                                                    gene_name + f'(exon{feature.qualifiers["number"][0]})',
                                                    feature.strand))
                    genes_in_location.append(gene_name)
            else:
                self.gene_locations.append(locations[i])
        self.gene_locations.append(locations[-1])  # 刚刚没有加入最后一个gene的位置
        self.gene_locations.sort(
            key=lambda location: location[0] - location[1] / (self.length + 1))  # 在起点相同的情况下，让终点靠后的基因排在前面

        # 计算向前平移的距离
        last_end = 0
        bias = 0
        for i in range(len(self.gene_locations)):
            if self.gene_locations[i][0] >= last_end:
                bias += int((self.gene_locations[i][0] - last_end) * (1 - interval_ratio))
            # else:
            #     assert False
            last_end = max(self.gene_locations[i][1], last_end)
            self.gene_location_bias.append(bias)


def my_circle_draw(inputfile, short_TR, long_TR, dis_R, outputfile, dis_R_Threshold, interval_ratio, draw_exon=True,
                   draw_shadow=True, draw_repeat_region=None, draw_misc_feature=None, draw_ncRNA=None,
                   draw_3UTR=None, draw_5UTR=None, draw_three_repeat=True, draw_ch_PAS=True,
                   draw_ch_PSB=True, draw_ch_PET=True, draw_ch_ATP=True, draw_ch_NDH=True, draw_ch_RBC=True,
                   draw_ch_RPO=True, draw_ch_RSP=True, draw_ch_RPL=True, draw_ch_CLP=True, draw_ch_YCF=True,
                   draw_ch_TRN=True, draw_ch_RRN=True, draw_ch_others=True):
    plt.rcParams['savefig.dpi'] = 3000  # 图片像素
    # plt.rcParams['figure.dpi'] = 300 # 分辨率
    feature_types_to_draw = ['CDS', 'tRNA', 'rRNA']

    # 读取 inputfile
    record = SeqIO.read(inputfile, "genbank")
    genes = genes_dict()

    for feature in record.features:
        genes.add_feature(feature)

    genes.compute_location_bias(interval_ratio)

    old_len = len(record.seq)
    new_len = int(old_len - genes.gene_locations[-1][1]) * interval_ratio + genes.gene_locations[-1][1] - \
              genes.gene_location_bias[-1]

    ratio = (2 * np.pi) / new_len
    ratio_old = (2 * np.pi) / old_len
    genes_hight = 14
    radius = 200
    label_font_size = 9
    label_height = 13
    font_size2height = {9: 13, 8: 10, 7: 10, 6: 9, 5: 6}  # label_font_size大小下文本的高度
    text_interval1 = 2 / 180 * np.pi  # 外圈文字角度差，2度
    text_interval2 = 3.5 / 180 * np.pi  # 内圈文字角度差，3.5度
    np.array(genes_hight)
    fig = plt.figure(figsize=(12, 12))  # 调节图象大小

    # plt.axis('off')
    # plt.gca().add_patch(
    #     plt.Rectangle(xy=(0, 0), width=0.1, height=0.1, fill=True, edgecolor='black', linewidth=0.5))
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6], polar=True)

    plt.axis('off')
    # x = plt.Circle((400, 400), 20, color='b', alpha=0.5)
    plt.polar(np.arange(0, 2 * np.pi, 2 * np.pi / 700), np.ones(700) * radius, color='black')  # 中轴线
    # plt.fill(np.arange(0, 2 * np.pi, 2 * np.pi / 700), np.ones(700) * 95, color=colors2, alpha=0.3)  # 背景色
    ax.xaxis.grid()
    ax.yaxis.grid()
    if 'organism' not in record.annotations.keys():
        title = f'species_name'
    else:
        title = record.annotations['organism']
    G_C_con = record.seq.count('G') + record.seq.count('C')
    lens = str('{:,}'.format(old_len)) + ' bp  GC: ' + str("%.2f" % (100 * G_C_con / len(record.seq))) + '%'

    if draw_three_repeat and interval_ratio == 1:
        plt.text(0.52, 0.95, title, ha='center', va='center', rotation=0,
                 wrap=True, fontsize=25, fontweight='bold', style='italic', transform=plt.gcf().transFigure)
        plt.text(0.52, 0.92, 'Plastid Genome', ha='center', va='center', rotation=0,
                 wrap=True, fontsize=17, transform=plt.gcf().transFigure)
        plt.text(0.52, 0.90, lens, ha='center', va='center', rotation=0,
                 wrap=True, fontsize=17, transform=plt.gcf().transFigure)
    else:
        plt.text(0.52, 0.53, title, ha='center', va='center', rotation=0,
                 wrap=True, fontsize=13, fontweight='bold', style='italic', transform=plt.gcf().transFigure)
        plt.text(0.52, 0.50, 'Plastid Genome', ha='center', va='center', rotation=0,
                 wrap=True, fontsize=13, transform=plt.gcf().transFigure)
        plt.text(0.52, 0.47, lens, ha='center', va='center', rotation=0,
                 wrap=True, fontsize=13, transform=plt.gcf().transFigure)

    dis_R = open(dis_R, 'r')  # 读取散在重复序列文件
    dis_location_old = dis_R.readlines()
    dis_location = []
    for i in range(len(dis_location_old)):
        dis_location.append((int(dis_location_old[i].replace('\n', '').split('\t')[0]),
                             int(dis_location_old[i].replace('\n', '').split('\t')[1]),
                             int(dis_location_old[i].replace('\n', '').split('\t')[2]),
                             int(dis_location_old[i].replace('\n', '').split('\t')[3])))

    if draw_shadow and interval_ratio == 1:
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
        if ratio_old * IR_A_start + 2 * np.pi - ratio_old * IR_B_end < ratio_old * IR_B_start - ratio_old * IR_A_end:
            theta_SSC, theta_LSC = theta_LSC, theta_SSC

        theta = np.linspace(theta_IR_A[0], theta_IR_A[1], 100)
        r = (radius + genes_hight) * np.ones_like(theta)
        plt.fill_between(theta, 0, r, color=(184 / 255, 211 / 255, 143 / 255), alpha=0.5)

        theta2 = np.linspace(theta_IR_B[0], theta_IR_B[1], 100)
        r2 = (radius + genes_hight) * np.ones_like(theta2)
        plt.fill_between(theta2, 0, r2, color=(184 / 255, 211 / 255, 143 / 255), alpha=0.5)

        theta3 = np.linspace(theta_LSC[0], theta_LSC[1], 100)
        r3 = (radius + genes_hight) * np.ones_like(theta3)
        plt.fill_between(theta3, 0, r3, color=colors[10], alpha=0.3)

        theta4 = np.linspace(theta_SSC[0], theta_SSC[1], 100)
        r4 = (radius + genes_hight) * np.ones_like(theta4)
        plt.fill_between(theta4, 0, r4, color=(238 / 255, 220 / 255, 130 / 255), alpha=0.4)

    thetas1 = []
    widths1 = []
    labels1 = []
    bar_colors1 = []
    thetas2 = []
    widths2 = []
    labels2 = []
    bar_colors2 = []
    thetas1_old = []
    thetas2_old = []
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
                            theta = ratio * ((start + end) / 2 - bias)
                            thetas1_old.append(theta)
                            thetas1.append(theta)
                            labels1.append(gene_name + f'(exon{feature.qualifiers["number"][0]})')
                            widths1.append(ratio * (end - start))
                            bar_colors1.append(colors[colors_index(gene_name)])
                        else:
                            theta = ratio * ((start + end) / 2 - bias)
                            thetas2_old.append(theta)
                            thetas2.append(theta)
                            labels2.append(gene_name + f'(exon{feature.qualifiers["number"][0]})')
                            widths2.append(ratio * (end - start))
                            bar_colors2.append(colors[colors_index(gene_name)])
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
                                    theta = ratio * ((start + end) / 2 - bias)
                                    thetas1_old.append(theta)
                                    thetas1.append(theta)
                                    labels1.append(gene_name)
                                    widths1.append(ratio * (end - start))
                                    bar_colors1.append(colors[colors_index(gene_name)])
                                else:
                                    theta = ratio * ((start + end) / 2 - bias)
                                    thetas2_old.append(theta)
                                    thetas2.append(theta)
                                    labels2.append(gene_name)
                                    widths2.append(ratio * (end - start))
                                    bar_colors2.append(colors[colors_index(gene_name)])
                    genes_flag[gene_name] = True
            else:
                if location[3] == 1:
                    theta = ratio * ((location[0] + location[1]) / 2 - genes.gene_location_bias[i])
                    thetas1_old.append(theta)
                    thetas1.append(theta)
                    labels1.append(gene_name)
                    widths1.append(ratio * (location[1] - location[0]))
                    bar_colors1.append(colors[colors_index(gene_name)])
                else:
                    theta = ratio * ((location[0] + location[1]) / 2 - genes.gene_location_bias[i])
                    thetas2_old.append(theta)
                    thetas2.append(theta)
                    labels2.append(gene_name)
                    widths2.append(ratio * (location[1] - location[0]))
                    bar_colors2.append(colors[colors_index(gene_name)])
    ax.bar(thetas1, genes_hight, width=widths1, bottom=radius, color=bar_colors1,
           edgecolor='black', linewidth=0.5)
    ax.bar(thetas2, -genes_hight, width=widths2, bottom=radius, color=bar_colors2,
           edgecolor='black', linewidth=0.5)

    renderer = ax.figure.canvas.get_renderer()

    # 排序
    indexs1 = np.arange(len(thetas1)).astype(int)
    for i in range(len(thetas1)):
        for j in range(len(thetas1) - 1):
            if thetas1[indexs1[j]] > thetas1[indexs1[j + 1]]:
                x = indexs1[j + 1]
                indexs1[j + 1] = indexs1[j]
                indexs1[j] = x
    thetas1 = [thetas1[i] for i in indexs1]
    thetas1_old = [thetas1_old[i] for i in indexs1]
    labels1 = [labels1[i] for i in indexs1]
    widths1 = [widths1[i] for i in indexs1]

    indexs2 = np.arange(len(thetas2))
    for i in range(len(thetas2)):
        for j in range(len(thetas2) - 1):
            if thetas2[indexs2[j]] > thetas2[indexs2[j + 1]]:
                x = indexs2[j + 1]
                indexs2[j + 1] = indexs2[j]
                indexs2[j] = x
    thetas2 = [thetas2[i] for i in indexs2]
    thetas2_old = [thetas2_old[i] for i in indexs2]
    labels2 = [labels2[i] for i in indexs2]
    widths2 = [widths2[i] for i in indexs2]

    # # 调整字体大小
    # font_indexs1 = label_font_size
    # indexs1_height = font_size2height[font_indexs1]
    # while len(indexs1) * indexs1_height / label_height * text_interval1 > 2 * np.pi and font_indexs1 > 1:
    #     font_indexs1 -= 1
    #     lab = ax.text(0, 0, 'A', transform=None, ha='center', va='center', fontsize=font_indexs1)
    #     bbox = lab.get_window_extent(renderer=renderer)
    #     lab.set_visible(False)
    #     indexs1_height = bbox.height
    # text_interval1 = indexs1_height / label_height * text_interval1
    #
    # font_indexs2 = label_font_size
    # indexs2_height = label_height
    # while len(indexs2) * indexs2_height / label_height * text_interval2 > 2 * np.pi and font_indexs2 > 1:
    #     font_indexs2 -= 1
    #     lab = ax.text(0, 0, 'A', transform=None, ha='center', va='center', fontsize=font_indexs2)
    #     bbox = lab.get_window_extent(renderer=renderer)
    #     lab.set_visible(False)
    #     indexs2_height = bbox.height
    # text_interval2 = indexs2_height / label_height * text_interval2

    # 去除重叠
    def remove1(thetas, font_size, standard_text_interval):
        thetas_origin = thetas.copy()
        text_interval = font_size2height[font_size] / label_height * standard_text_interval
        thetas_old = thetas.copy()
        index2cluster_old = np.arange(len(thetas), dtype=int)
        index2cluster_new = index2cluster_old.copy()
        while True:
            for i in range(len(thetas) - 1):
                if thetas[i + 1] - thetas[i] < text_interval - 1e-6:
                    new_cluster = index2cluster_new[i]
                    old_cluster = index2cluster_old[i + 1]
                    if new_cluster != old_cluster and font_size > 5 and 1 + (
                            index2cluster_new == new_cluster).sum() > 30:
                        return remove1(thetas_origin, font_size - 1, standard_text_interval)
                    else:
                        # index2cluster_new[index2cluster_old == old_cluster] = new_cluster
                        index2cluster_new[i + 1] = index2cluster_new[i]
                        thetas[i + 1] = thetas[i] + text_interval
                        center = (thetas[int(np.where(index2cluster_new == new_cluster)[0][0])] +
                                  thetas[int(np.where(index2cluster_new == new_cluster)[0][-1])]) / 2

                        cluster_width = len(np.where(index2cluster_new == new_cluster)[0]) * text_interval
                        target_center = min(np.pi * 2 - cluster_width / 2, max(cluster_width / 2, (
                                thetas_origin[int(np.where(index2cluster_new == new_cluster)[0][0])] +
                                thetas_origin[int(np.where(index2cluster_new == new_cluster)[0][-1])]) / 2))
                        for j in np.where(index2cluster_new == new_cluster)[0]:
                            thetas[j] += target_center - center

            if (index2cluster_new == index2cluster_old).all():
                break

            # clusters = []
            # for i in range(len(thetas)):
            #     cluster = index2cluster_new[i]
            #     if cluster not in clusters:
            #         clusters.append(cluster)
            #         center = (thetas[int(np.where(index2cluster_new == cluster)[0][0])] +
            #                   thetas[int(np.where(index2cluster_new == cluster)[0][-1])]) / 2
            #
            #         cluster_width = len(np.where(index2cluster_new == cluster)[0]) * text_interval
            #         target_center = min(np.pi * 2 - cluster_width / 2, max(cluster_width / 2, (
            #                 thetas_old[int(np.where(index2cluster_new == cluster)[0][0])] +
            #                 thetas_old[int(np.where(index2cluster_new == cluster)[0][-1])]) / 2))
            #         for j in np.where(index2cluster_new == cluster)[0]:
            #             thetas[j] += target_center - center
            index2cluster_old = index2cluster_new.copy()
        return thetas, font_size

    thetas1, font_indexs1 = remove1(thetas1, label_font_size, text_interval1)

    # 加外圈标签
    for i in range(len(thetas1)):
        theta_old, theta, width, label = thetas1_old[i], thetas1[i], widths1[i], labels1[i]
        rotation = np.rad2deg(theta)
        lab = ax.text(0, 0, label, transform=None, ha='center', va='center', fontsize=font_indexs1)
        bbox = lab.get_window_extent(renderer=renderer)
        ax.plot([theta_old, theta_old], [214.5, 217], color="black", linewidth=1)
        ax.plot([theta_old, theta], [217, 222], color="black", linewidth=1)
        ax.plot([theta, theta], [222, 224.5], color="black", linewidth=1)
        if np.pi / 2 < theta < 3 * np.pi / 2:
            lab.set_position((theta, .97 + bbox.width / (2 * 370)))
            lab.set_transform(ax.get_xaxis_transform())
            lab.set_rotation(rotation + 180)
        else:
            lab.set_position((theta, .97 + bbox.width / (2 * 370)))
            lab.set_transform(ax.get_xaxis_transform())
            lab.set_rotation(rotation)

    thetas2, font_indexs2 = remove1(thetas2, label_font_size, text_interval2)
    # 加内圈标签
    for i in range(len(thetas2)):
        # for bar, theta, width, strand, label in zip(bars, thetas, widths, strands, labels):

        theta, theta_old, width, label = thetas2[i], thetas2_old[i], widths2[i], labels2[i]
        rotation = np.rad2deg(theta)
        lab = ax.text(0, 0, label, transform=None, ha='center', va='center', fontsize=font_indexs2)
        bbox = lab.get_window_extent(renderer=renderer)
        # offset = (180 + bars[0].get_height()-bbox.width)/(y1-y0) # 调节标签离心位置
        # invb = ax.transData.inverted().transform([[0, 0],[bbox.width, 0] ])
        ax.plot([theta_old, theta_old], [185.8, 183.3], color="black", linewidth=1)
        ax.plot([theta_old, theta], [183.3, 178.3], color="black", linewidth=1)
        ax.plot([theta, theta], [178.3, 175.8], color="black", linewidth=1)
        if np.pi / 2 < theta < 3 * np.pi / 2:
            lab.set_position((theta, 0.75 - bbox.width / (2 * 340)))
            lab.set_transform(ax.get_xaxis_transform())
            lab.set_rotation(rotation + 180)
        else:
            lab.set_position((theta, 0.75 - bbox.width / (2 * 340)))
            lab.set_transform(ax.get_xaxis_transform())
            lab.set_rotation(rotation)

    if draw_three_repeat and interval_ratio == 1:
        plt.polar(np.arange(0, 2 * np.pi, 2 * np.pi / 2000), np.ones(2000) * 93, color='black', lw=1)
        plt.polar(np.arange(0, 2 * np.pi, 2 * np.pi / 700), np.ones(700) * 89, color='black', lw=1)  # GC含量中间虚线
        plt.polar(np.arange(0, 2 * np.pi, 2 * np.pi / 700), np.ones(700) * 85, ls='--', color='black', lw=1)  # GC含量线
        plt.polar(np.arange(0, 2 * np.pi, 2 * np.pi / 700), np.ones(700) * 77, color='black', lw=1)  # 微卫星重复实线
        plt.polar(np.arange(0, 2 * np.pi, 2 * np.pi / 700), np.ones(700) * 69, color='black', lw=1)  # 串联重复实线
        plt.polar(np.arange(0, 2 * np.pi, 2 * np.pi / 700), np.ones(700) * 61, color='black', lw=1)  # 串联重复实线
        # GC含量
        GC_content = []
        for i in range(int(old_len / 100)):
            G_content = record.seq[100 * i:(i + 1) * 100].count('G')
            C_content = record.seq[100 * i:(i + 1) * 100].count('C')
            GC_content.append(G_content + C_content)

        # 画GC含量直方图
        width = 2 * np.pi / len(GC_content)
        for i in range(len(GC_content)):
            ax.plot([i * width, i * width], [90, (GC_content[i] - 50) * 10 / 100 + 90], color="orange", linewidth=0.5)

        # 画刻度
        for i in range(12):
            theta = old_len / 12 * i * ratio_old
            rotation = np.rad2deg(theta)
            scale_label = str(int(old_len / 12000 * i)) + 'kb'
            lab = ax.text(0, 0, scale_label, transform=None, ha='center', va='center', color='blue',fontsize=8)
            bbox = lab.get_window_extent(renderer=renderer)
            ax.plot([theta, theta], [93, 95], color="black", linewidth=1)
            if 0 < theta < np.pi:
                lab.set_position((theta, (93 - bbox.height / 2 * 0.95) / 260 + 0.078))
                lab.set_transform(ax.get_xaxis_transform())
                lab.set_rotation(rotation + 270)
            else:
                lab.set_position((theta, (93 - bbox.height / 2 * 0.95) / 260 + 0.078))
                lab.set_transform(ax.get_xaxis_transform())
                lab.set_rotation(rotation + 90)

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

        for i in short_location:
            theta = ratio_old * (i[0] + i[1]) / 2
            # ax.bar(theta, 14, width=ratio_old * (i[1] - i[0]), bottom=193, color='Tomato', edgecolor='Tomato', linewidth=1.5, alpha=0.6)
            ax.bar(theta, 7.5, width=ratio_old * (i[1] - i[0]), bottom=77.5, color='Tomato', edgecolor='Tomato',
                   linewidth=1.5)

        for i in long_location:
            theta = ratio_old * (i[0] + i[1]) / 2
            # plt.scatter(theta, 200, s=40, c='BlueViolet', marker='v', edgecolor='purple', alpha=0.8)
            # ax.bar(theta, 14, width=ratio_old * (i[1] - i[0]), bottom=193, color='BlueViolet', edgecolor='BlueViolet', linewidth=1.5, alpha=0.6)
            ax.bar(theta, 7.5, width=ratio_old * (i[1] - i[0]), bottom=69.5, color='BlueViolet', edgecolor='BlueViolet',
                   linewidth=1.5)

        r1 = 15
        r2 = 60
        r3 = 50
        r4 = 200
        r5 = 70
        for i in range(len(dis_location)):
            thetas3 = []
            r3_s = []
            theta1 = ratio_old * (dis_location[i][1] + dis_location[i][0] / 2)
            theta2 = ratio_old * (dis_location[i][3] + dis_location[i][0] / 2)
            if theta1 > theta2:
                theta2, theta1 = theta1, theta2
            if theta1 - theta2 > np.pi:
                theta2 += 2 * np.pi
            elif theta2 - theta1 > np.pi:
                theta1 += 2 * np.pi
            if theta1 > theta2:
                theta1, theta2 = theta2, theta1
            theta4 = theta1 - ratio_old * dis_location[i][0] / 2
            theta5 = theta1 + ratio_old * dis_location[i][0] / 2
            theta6 = theta2 - ratio_old * dis_location[i][0] / 2
            theta7 = theta2 + ratio_old * dis_location[i][0] / 2
            if theta2 - theta1 < np.pi:
                # theta4 = theta1 - ratio_old * dis_location[i][0] / 2
                # theta5 = theta1 + ratio_old * dis_location[i][0] / 2
                # theta6 = theta2 - ratio_old * dis_location[i][0] / 2
                # theta7 = theta2 + ratio_old * dis_location[i][0] / 2
                bias_47 = r3
            else:
                # theta4 = theta1 + ratio_old * dis_location[i][0] / 2
                # theta5 = theta1 - ratio_old * dis_location[i][0] / 2
                # theta6 = theta2 + ratio_old * dis_location[i][0] / 2
                # theta7 = theta2 - ratio_old * dis_location[i][0] / 2
                bias_47 = r5
            # theta4、theta5为第一个重复序列的两端， theta6、theta7位第二个重复序列的两端
            if dis_location[i][0] > dis_R_Threshold:
                theta_line = []
                r_line = []

                # theta4到theta5
                theta_line.append(theta4)
                r_line.append(radius)
                while theta_line[-1] < theta5 - 0.05:
                    theta_line.append(theta_line[-1] + 0.05)
                    r_line.append(radius)
                theta_line.append(theta5)
                r_line.append(radius)

                # theta5到theta6
                if r4 * np.cos((theta5 - theta6) / 2) < r3:
                    theta_line.append(theta6)
                    r_line.append(radius)
                else:
                    theta_line1, r_line1 = draw_parabola(r5 + r3 - bias_47, r4, theta5, theta6)
                    theta_line += theta_line1
                    r_line += r_line1

                # theta6到theta7:
                theta_line.append(theta6)
                r_line.append(radius)
                while theta_line[-1] < theta7 - 0.05:
                    theta_line.append(theta_line[-1] + 0.05)
                    r_line.append(radius)
                theta_line.append(theta7)
                r_line.append(radius)

                # theta7到theta4
                if r4 * np.cos((theta7 - theta4) / 2) < r3:
                    theta_line.append(theta4)
                    r_line.append(radius)
                else:
                    theta_line1, r_line1 = draw_parabola(bias_47, r4, theta4, theta7)
                    theta_line += theta_line1[::-1]
                    r_line += r_line1

                # if dis_location[i][2] == 1:
                #     ax.bar(theta1, 10, width=ratio_old * dis_location[i][0], bottom=195, color='LightSalmon',
                #            edgecolor='grey', linewidth=1.5, alpha=0.5)
                #     ax.bar(theta2, 10, width=ratio_old * dis_location[i][0], bottom=195, color='LightSalmon',
                #            edgecolor='grey', linewidth=1.5, alpha=0.5)
                #
                #     plt.fill(theta_line, r_line, color='LightSalmon', alpha=0.3,
                #              edgecolor='grey')
                #     # if r4 * np.cos((theta2 - theta1) / 2) < r3:
                #     #     a = [theta4, (theta4 + theta5) / 2, theta5, theta6, (theta6 + theta7) / 2, theta7]
                #     #     b = [radius, radius, radius, radius, radius, radius]
                #     #     plt.fill(a, b, color='LightSalmon', alpha=0.3, edgecolor='grey')  # 背景色
                #     # else:
                #     #     theta_line1, r_line1 = draw_parabola(r5, r4, theta5, theta6)
                #     #     theta_line2, r_line2 = draw_parabola(r3, r4, theta7, theta4)
                #     #     plt.fill(theta_line1 + theta_line2[::-1], r_line1 + r_line2[::-1], color='LightSalmon', alpha=0.3,
                #     #              edgecolor='grey')  # 背景色
                # if dis_location[i][2] == -1:
                #     ax.bar(theta1, 10, width=ratio_old * dis_location[i][0], bottom=195, color='Aquamarine',
                #            edgecolor='grey', linewidth=1.5, alpha=0.5)
                #     ax.bar(theta2, 10, width=ratio_old * dis_location[i][0], bottom=195, color='Aquamarine',
                #            edgecolor='grey', linewidth=1.5, alpha=0.5)
                #     # if r4 * np.cos((theta2 - theta1) / 2) < r3:
                #     #     a = [theta4, (theta4 + theta5) / 2, theta5, theta6, (theta6 + theta7) / 2, theta7]
                #     #     b = [radius, radius, radius, radius, radius, radius]
                #     #     plt.fill(a, b, color='Aquamarine', alpha=0.3, edgecolor='grey')  # 背景色
                #     # else:
                #     #     theta_line1, r_line1 = draw_parabola(r5, r4, theta5, theta6)
                #     #     theta_line2, r_line2 = draw_parabola(r3, r4, theta7, theta4)
                #     plt.fill(theta_line, r_line, color='Aquamarine', alpha=0.3,
                #                  edgecolor='grey')  # 背景色
            if dis_location[i][2] == 1:
                ax.bar(theta1, 7.5, width=ratio_old * dis_location[i][0], bottom=61.5, color='orange',
                       edgecolor='orange',
                       linewidth=1.5, alpha=0.6)
                ax.bar(theta2, 7.5, width=ratio_old * dis_location[i][0], bottom=61.5, color='orange',
                       edgecolor='orange',
                       linewidth=1.5, alpha=0.6)
            if dis_location[i][2] == -1:
                ax.bar(theta1, 7.5, width=ratio_old * dis_location[i][0], bottom=61.5, color='green', edgecolor='green',
                       linewidth=1.5, alpha=0.6)
                ax.bar(theta2, 7.5, width=ratio_old * dis_location[i][0], bottom=61.5, color='green', edgecolor='green',
                       linewidth=1.5, alpha=0.6)

            if theta2 == theta1:
                print(1)
            if r2 * np.cos((theta2 - theta1) / 2) < r1:
                ax.plot([theta1, theta2], [r2, r2])
            else:
                ax.plot(draw_parabola(r1, r2, theta1, theta2)[0], draw_parabola(r1, r2, theta1, theta2)[1], alpha=0.7)

    for i in range(1, len(record.features)):
        if draw_repeat_region:
            if record.features[i].type == 'repeat_region':
                theta = ratio_old * (int(record.features[i].location.start) + int(record.features[i].location.end)) / 2
                plt.scatter(theta, 194, s=40, edgecolors=(102 / 255, 205 / 255, 170 / 255), marker='v', alpha=0.8,
                            linewidths=1.3, facecolors='none')
        if draw_misc_feature:
            if record.features[i].type == 'misc_feature':
                theta = ratio_old * (int(record.features[i].location.start) + int(record.features[i].location.end)) / 2
                plt.scatter(theta, 200, s=40, edgecolors=(226 / 255, 0 / 255, 26 / 255), marker='o', alpha=0.8, linewidths=1.3,
                            facecolors='none')
        if draw_ncRNA:
            if record.features[i].type == 'ncRNA':
                theta = ratio_old * (int(record.features[i].location.start) + int(record.features[i].location.end)) / 2
                plt.scatter(theta, 206, s=40, edgecolors=(255 / 255, 105 / 255, 180 / 255), marker='*', alpha=0.8,
                            linewidths=1.3, facecolors='none')
        if draw_5UTR:
            if record.features[i].type == "5'UTR":
                theta = ratio_old * (int(record.features[i].location.start) + int(record.features[i].location.end)) / 2
                plt.scatter(theta, 188, s=40, edgecolors=(151 / 255, 190 / 255, 13 / 255), marker='P', alpha=0.8, linewidths=1.3,
                            facecolors='none')
        if draw_3UTR:
            if record.features[i].type == "3'UTR":
                theta = ratio_old * (int(record.features[i].location.start) + int(record.features[i].location.end)) / 2
                plt.scatter(theta, 212, s=40, edgecolors='brown', marker='d', alpha=0.8, linewidths=1.3,
                            facecolors='none')

    ax.set_rlim(bottom=0, top=235)
    # ax.set_rmax(np.max(224.5))
    # 添加图例 左边
    ax = fig.add_axes([0.2, -0.02, 0.5, 0.11],
                      polar=False)  # 增加一个新画布（注意调整画布大小和范围，避免覆盖之前的极坐标画布），注意polar要选False，选了True的话这个画布的坐标就是极坐标
    plt.axis('off')

    plt.gca().add_patch(
        plt.Rectangle(xy=(0, 0.62), width=0.014, height=0.05, fill=True, color=colors[0]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.001, 0.62), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.03, 0.64, 'photosystem I', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0, 0.54), width=0.014, height=0.05, fill=True, color=colors[1]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.001, 0.54), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.03, 0.56, 'photosystem II', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0, 0.46), width=0.014, height=0.05, fill=True, color=colors[2]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.001, 0.46), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.03, 0.48, 'cytochrome b/f complex', ha='left', va='center', rotation=0,
             wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0, 0.38), width=0.014, height=0.05, fill=True, color=colors[3]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.001, 0.38), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.03, 0.40, 'ATP synthesis', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0, 0.30), width=0.014, height=0.05, fill=True, color=colors[4]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.001, 0.30), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.03, 0.32, 'NADH dehydrogenase', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0, 0.22), width=0.014, height=0.05, fill=True, color=colors[5]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.001, 0.22), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.03, 0.24, 'Rubisco large subunit', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0, 0.14), width=0.014, height=0.05, fill=True, color=colors[6]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.001, 0.14), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.03, 0.16, 'RNA polymerase', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0, 0.06), width=0.014, height=0.05, fill=True, color=colors[7]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.001, 0.06), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.03, 0.08, 'small ribosomal protein', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0, -0.02), width=0.014, height=0.05, fill=True, color=colors[8]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.001, -0.02), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.03, 0, 'large ribosomal protein', ha='left', va='center', rotation=0, wrap=True)

    plt.text(0, -0.08, 'The colored parabola in the center circle represents the Dispersed Repeats', ha='left',
             va='center', rotation=0, wrap=True)

    # 添加图例 右边
    plt.gca().add_patch(
        plt.Rectangle(xy=(0.33, 0.62), width=0.014, height=0.05, fill=True, color=colors[9]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.331, 0.62), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.36, 0.64, 'clpP,matK,infA', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0.33, 0.54), width=0.014, height=0.05, fill=True, color=colors[10]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.331, 0.54), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.36, 0.56, 'hypothetical reading frame', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0.33, 0.46), width=0.014, height=0.05, fill=True, color=colors[11]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.331, 0.46), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.36, 0.48, 'ribosomal RNA', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0.33, 0.38), width=0.014, height=0.05, fill=True, color=colors[12]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.331, 0.38), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.36, 0.40, 'transfer RNA', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0.33, 0.30), width=0.014, height=0.05, fill=True, color=colors[13]))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.331, 0.30), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.36, 0.32, 'other', ha='left', va='center', rotation=0, wrap=True)

    plt.gca().add_patch(
        plt.Rectangle(xy=(0.33, 0.22), width=0.014, height=0.05, fill=True, color='Tomato'))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.331, 0.22), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.36, 0.24, 'Short Tandem Repeats', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0.33, 0.14), width=0.014, height=0.05, fill=True, color='BlueViolet'))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.331, 0.14), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.36, 0.16, 'Long Tandem Repeats', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0.33, 0.06), width=0.014, height=0.05, fill=True, color='orange', alpha=0.6))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.331, 0.06), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.36, 0.08, 'Dispersed Repeats(direct matches)', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(
        plt.Rectangle(xy=(0.33, -0.02), width=0.014, height=0.05, fill=True, color='green', alpha=0.6))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.331, -0.02), width=0.014, height=0.05, fill=False, edgecolor='black',
                                      linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.36, 0, 'Dispersed Repeats(palindromic matches)', ha='left', va='center', rotation=0, wrap=True)



    plt.scatter(0.75, 0.64, s=60, edgecolors=(226 / 255, 0 / 255, 26 / 255), marker='o', alpha=0.8, linewidths=1.3, facecolors='none')
    plt.text(0.77, 0.64, 'misc feature', ha='left', va='center', rotation=0, wrap=True)
    plt.scatter(0.75, 0.56, s=60, edgecolors=(255 / 255, 105 / 255, 180 / 255), marker='*', alpha=0.8, linewidths=1.3,
                facecolors='none')
    plt.text(0.77, 0.56, 'ncRNA', ha='left', va='center', rotation=0, wrap=True)
    plt.scatter(0.75, 0.48, s=60, edgecolors=(102 / 255, 205 / 255, 170 / 255), marker='v', alpha=0.8, linewidths=1.3,
                facecolors='none')
    plt.text(0.77, 0.48, 'repeat region', ha='left', va='center', rotation=0, wrap=True)
    plt.scatter(0.75, 0.40, s=60, edgecolors=(151 / 255, 190 / 255, 13 / 255), marker='P', alpha=0.8, linewidths=1.3, facecolors='none')
    plt.text(0.77, 0.40, "5'UTR", ha='left', va='center', rotation=0, wrap=True)
    plt.scatter(0.75, 0.32, s=60, edgecolors='brown', marker='d', alpha=0.8, linewidths=1.3, facecolors='none')
    plt.text(0.77, 0.32, "3'UTR", ha='left', va='center', rotation=0, wrap=True)

    plt.gca().add_patch(plt.Rectangle(xy=(0.745,0.22), width=0.014,height=0.05, fill=True, color=(184/255,211/255,143/255), alpha=0.6))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.746, 0.22), width=0.014, height=0.05, fill=False, edgecolor='black',linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.775, 0.24, 'Inverted Repeat(IR)', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(plt.Rectangle(xy=(0.745,0.14), width=0.014,height=0.05, fill=True, color=colors[10], alpha=0.4))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.746, 0.14), width=0.014, height=0.05, fill=False, edgecolor='black',linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.775, 0.16, 'Large Single Copy(LSR)', ha='left', va='center', rotation=0, wrap=True)
    plt.gca().add_patch(plt.Rectangle(xy=(0.745,0.06), width=0.014,height=0.05, fill=True, color=(238/255, 220/255, 130/255), alpha=0.5))  # xy是矩形的左下角坐标
    plt.gca().add_patch(plt.Rectangle(xy=(0.746, 0.06), width=0.014, height=0.05, fill=False, edgecolor='black',linewidth=0.5))  # xy是矩形的左下角坐标
    plt.text(0.775, 0.08, 'Small Single Copy(SSR)', ha='left', va='center', rotation=0, wrap=True)
    plt.savefig(outputfile, bbox_inches='tight', pad_inches=0.5, dpi=200, transparent=True)
    #plt.show()


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
    parser.add_argument('--interval_ratio', help='the ratio of interval to be kept,0-1', type=float, default=0.1)
    parser.add_argument('--draw_exon', help='the ratio of interval to be kept,0-1', type=bool, default=True)
    parser.add_argument('--draw_shadow', help='the ratio of interval to be kept,0-1', type=bool, default=True)
    parser.add_argument('--draw_ch_NDH', help='the ratio of interval to be kept,0-1', type=bool, default=False)
    parser.add_argument('--draw_ch_PAS', help='the ratio of interval to be kept,0-1', type=bool, default=False)
    parser.add_argument('--draw_ch_ATP', help='the ratio of interval to be kept,0-1', type=bool, default=False)
    args = parser.parse_args()
    my_circle_draw(args.input, args.short_TR, args.long_TR, args.dispersed_R, args.output, args.dis_R_Threshold,
                   args.interval_ratio, args.draw_exon, args.draw_shadow,
                   draw_ch_NDH=args.draw_ch_NDH,draw_ch_PAS=args.draw_ch_PAS,draw_ch_ATP=args.draw_ch_ATP, draw_ch_others=False)
