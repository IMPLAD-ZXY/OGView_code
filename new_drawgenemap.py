import argparse

import numpy as np
from Bio import SeqIO
from PIL import Image, ImageDraw, ImageFont


def mydrawgenemap(inputfile, outputfile):
    record = SeqIO.read(inputfile, "genbank")  # 读取gb文件

    colors = [(241, 204, 184), (241, 241, 184), (184, 211, 143), (184, 241, 237), (217, 184, 241),
              (221, 255, 149), (184, 241, 204), (231, 219, 202), (207, 136, 120), (239, 87, 103), (196, 144, 160),
              (255, 150, 113), (255, 199, 95), (249, 248, 113),
              (221, 255, 149), (184, 241, 204), (132, 94, 194), (214, 93, 177), (255, 111, 145), (231, 219, 202),
              (207, 136, 120), (239, 87, 103), (196, 144, 160), (132, 94, 194), (214, 93, 177), (255, 111, 145),
              (255, 150, 113), (255, 199, 95), (249, 248, 113),
              (221, 255, 149), (184, 241, 204), (231, 219, 202), (207, 136, 120), (239, 87, 103), (196, 144, 160)]
    color2 = (22, 41, 131)
    total_width = 5000
    width_ratio = total_width / 10000
    total_height = int(5000 * width_ratio)
    blank = int(300 * width_ratio)  # 中央横线到两边的距离
    text_blank = int(300 * width_ratio)  # 文本到两边的最大距离
    ratio = (total_width - blank * 2) / len(record.seq)
    img = Image.new('RGB', (total_width, total_height), 'white')
    hight = int(200 * width_ratio)  # 每个矩形色块的高度
    line_width = int(40 * width_ratio)  # 中轴线粗度
    line_with2 = int(7.5 * width_ratio)  # 延伸线的粗细
    font_size = int(150 * width_ratio)
    font_style = '/home/xyzhang/pmgdraw/bin/arial.ttf'
    font = ImageFont.truetype(font_style, font_size)
    y_re1 = (total_height - line_width) / 2  # 横线上方的色块
    y_re2 = (total_height + line_width) / 2  # 横线下方的色块
    y_text1 = (total_height - line_width) / 2 - hight - font_size - int(300 * width_ratio)  # 横线上方的exon字体
    y_text2 = (total_height + line_width) / 2 + hight + int(420 * width_ratio)  # 横线下方的exon字体
    text_interval = int(100 * width_ratio)
    imgDraw = ImageDraw.Draw(img)
    imgDraw.line((blank, total_height / 2, total_width - blank, total_height / 2), fill='black',
                 width=line_width)  # 中央横线

    genes = {}
    for i in range(len(record.features)):
        if record.features[i].type[:4] == 'exon':
            gene_name = record.features[i].qualifiers['gene'][0]
            if gene_name not in genes.keys():
                genes[gene_name] = record.features[i].location

    coords = {}
    for i in range(len(record.features)):
        if record.features[i].type == 'coord':
            gene_name = record.features[i].qualifiers['gene'][0]
            if gene_name not in coords.keys():
                coords[gene_name] = record.features[i].location

    # 调整字体大小
    nag_text = []
    pos_text = []
    for i in genes.keys():
        if genes[i].strand == 1:
            pos_text.append(i)
        else:
            nag_text.append(i)
    while True:
        nag_len = -text_interval
        pos_len = -text_interval
        for i in nag_text:
            nag_len += font.getsize(i)[0] + text_interval

        for i in pos_text:
            pos_len += font.getsize(i)[0] + text_interval
        if nag_len <= (total_width - blank * 2) and pos_len <= (total_width - blank * 2):
            break
        else:
            font_size -= 5
            font = ImageFont.truetype(font_style, font_size)

    label_location = []
    for i in genes.keys():
        if genes[i].strand == 1:
            if '(' in i:
                color = colors[int(i[4:i.index('(')]) - 1]
            else:
                color = colors[int(i[4:]) - 1]
            imgDraw.rectangle(
                (blank + ratio * genes[i].nofuzzy_start, y_re1 - hight, blank + ratio * genes[i].nofuzzy_end, y_re1),
                fill=color,
                outline='black', width=int(5 * width_ratio))  # 右上角，左下角 #横线上方的色块
            # imgDraw.text((blank+ratio*(genes[i].nofuzzy_start+genes[i].nofuzzy_end)/2, y_text1), text=i, fill='black',
            #      font=ImageFont.truetype('D:/Users/Administrator/PycharmProjects/pythonProject/新建文件夹/arial.ttf', font),
            #      spacing=5, anchor='mm')# 横线上方的exon
            left = blank + ratio * (genes[i].nofuzzy_start + genes[i].nofuzzy_end) / 2 - font.getsize(i)[0] / 2
            right = blank + ratio * (genes[i].nofuzzy_start + genes[i].nofuzzy_end) / 2 + font.getsize(i)[0] / 2
            label_location.append([left, right, (left + right) / 2, i, 1])
        elif genes[i].strand == -1:
            if '(' in i:
                color = colors[int(i[4:i.index('(')]) - 1]
            else:
                color = colors[int(i[4:]) - 1]
            imgDraw.rectangle(
                (blank + ratio * genes[i].nofuzzy_start, y_re2 + hight, blank + ratio * genes[i].nofuzzy_end, y_re2),
                fill=color,
                outline='black', width=int(5 * width_ratio))  # 横线下方的色块
            # imgDraw.text((blank+ratio*(genes[i].nofuzzy_start+genes[i].nofuzzy_end)/2, y_text2), text=i, fill='black',
            #          font=ImageFont.truetype('D:/Users/Administrator/PycharmProjects/pythonProject/新建文件夹/arial.ttf', font),
            #          spacing=5, anchor='mm')# 横线下方的exon
            left = blank + ratio * (genes[i].nofuzzy_start + genes[i].nofuzzy_end) / 2 - font.getsize(i)[0] / 2
            right = blank + ratio * (genes[i].nofuzzy_start + genes[i].nofuzzy_end) / 2 + font.getsize(i)[0] / 2
            label_location.append([left, right, (left + right) / 2, i, -1])

    for i in coords.keys():
        if coords[i].strand == 1:
            imgDraw.line((blank + ratio * coords[i].nofuzzy_start, y_re1 - hight,
                          blank + ratio * coords[i].nofuzzy_start, y_re1), fill=color2,
                         width=int(10 * width_ratio))  # 右上角，左下角 #横线上方的色块
            # imgDraw.text((blank+ratio*(coords[i].nofuzzy_start+coords[i].nofuzzy_end)/2, y_text1), text=i, fill='black',
            #      font=ImageFont.truetype('D:/Users/Administrator/PycharmProjects/pythonProject/新建文件夹/arial.ttf', font),
            #      spacing=5, anchor='mm')# 横线上方的exon
            left = blank + ratio * (coords[i].nofuzzy_start + coords[i].nofuzzy_end) / 2 - font.getsize(i)[0] / 2
            right = blank + ratio * (coords[i].nofuzzy_start + coords[i].nofuzzy_end) / 2 + font.getsize(i)[0] / 2
            label_location.append([left, right, (left + right) / 2, i, 1])  # 左、右、对应的居中位置、文本。注意居中位置在后续中不改变，是原先的
        elif coords[i].strand == -1:
            imgDraw.line((blank + ratio * coords[i].nofuzzy_start, y_re2 + hight,
                          blank + ratio * coords[i].nofuzzy_start, y_re2), fill=color2,
                         width=int(10 * width_ratio))  # 右上角，左下角 #横线下方的色块
            # imgDraw.text((blank+ratio*(coords[i].nofuzzy_start+coords[i].nofuzzy_end)/2, y_text2), text=i, fill='black',
            #          font=ImageFont.truetype('D:/Users/Administrator/PycharmProjects/pythonProject/新建文件夹/arial.ttf', font),
            #          spacing=5, anchor='mm')# 横线下方的exon
            left = blank + ratio * (coords[i].nofuzzy_start + coords[i].nofuzzy_end) / 2 - font.getsize(i)[0] / 2
            right = blank + ratio * (coords[i].nofuzzy_start + coords[i].nofuzzy_end) / 2 + font.getsize(i)[0] / 2
            label_location.append([left, right, (left + right) / 2, i, -1])  # 左、右、对应的居中位置、文本。注意居中位置在后续中不改变，是原先的

    # 解决重叠

    # 排序

    for i in range(len(label_location)):
        for index in range(len(label_location) - 1):
            if label_location[index][2] > label_location[index + 1][2]:
                value = label_location.pop(index)
                label_location.insert(index + 1, value)

    # for i in range(len(label_location)-1):
    #     if label_location[i+1][0] - label_location[i][1] < 10:
    #         overlap.append(label_location[i])
    #         overlap.append(label_location[i+1])
    #
    # while i < len(label_location):
    #     for j in overlap:
    #         if label_location[i] == j:
    #             label_location.pop(i)

    def remove_overlap(data):
        index2cluster_old = np.arange(len(data), dtype=int)
        while True:
            index2cluster_new = index2cluster_old.copy()
            for i in range(len(data) - 1):
                if data[i + 1][0] - data[i][1] < text_interval:
                    new_cluster = index2cluster_new[i]
                    old_cluster = index2cluster_old[i + 1]
                    index2cluster_new[index2cluster_old == old_cluster] = new_cluster
                    data[i + 1][1] = data[i + 1][1] - data[i + 1][0] + data[i][1] + text_interval
                    data[i + 1][0] = data[i][1] + text_interval

            if (index2cluster_new == index2cluster_old).all():
                break

            clusters = []
            for i in range(len(data)):
                cluster = index2cluster_new[i]
                if cluster not in clusters:
                    clusters.append(cluster)
                    center = (data[int(np.where(index2cluster_new == cluster)[0][0])][0] +
                              data[int(np.where(index2cluster_new == cluster)[0][-1])][1]) / 2

                    cluster_width = data[int(np.where(index2cluster_new == cluster)[0][-1])][1] - \
                                    data[int(np.where(index2cluster_new == cluster)[0][0])][0]
                    target_center = min(total_width - text_blank - cluster_width / 2,
                                        max(text_blank + cluster_width / 2,
                                            (data[int(np.where(index2cluster_new == cluster)[0][0])][2] +
                                             data[int(np.where(index2cluster_new == cluster)[0][-1])][2]) / 2))
                    for j in np.where(index2cluster_new == cluster)[0]:
                        data[j][0] += target_center - center
                        data[j][1] += target_center - center
            index2cluster_old = index2cluster_new.copy()
        return data

    label_location = remove_overlap(label_location)
    for i in range(len(label_location)):
        left, right, center, text, strand = label_location[i]
        if strand == 1:
            imgDraw.text(((left + right) / 2, y_text1), text=text, fill='black',
                         font=font, spacing=5, anchor='mm')  # 横线上方的文字标注
            imgDraw.line((center, y_re1 - hight, center, y_re1 - hight - int(50 * width_ratio)), fill='grey',
                         width=line_with2)  # line1
            imgDraw.line(
                (center, y_re1 - hight - int(50 * width_ratio), (left + right) / 2, y_text1 + int(150 * width_ratio)),
                fill='grey',
                width=line_with2)  # line2
            imgDraw.line(((left + right) / 2, y_text1 + int(150 * width_ratio), (left + right) / 2,
                          y_text1 + int(100 * width_ratio)), fill='grey',
                         width=line_with2)  # line3

        else:
            imgDraw.text(((left + right) / 2, y_text2), text=text, fill='black',
                         font=font, spacing=5, anchor='mm')  # 横线下方的文字标注
            imgDraw.line((center, y_re2 + hight, center, y_re2 + hight + int(50 * width_ratio)), fill='grey',
                         width=line_with2)  # line1
            imgDraw.line(
                (center, y_re2 + hight + int(50 * width_ratio), (left + right) / 2, y_text2 - int(150 * width_ratio)),
                fill='grey',
                width=line_with2)  # line2
            imgDraw.line(((left + right) / 2, y_text2 - int(150 * width_ratio), (left + right) / 2,
                          y_text2 - int(100 * width_ratio)), fill='grey',
                         width=line_with2)  # line3
    img.save(outputfile)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='the RNA file to edit', type=str,
                        default='D:/Users/Administrator/PycharmProjects/pythonProject/新建文件夹/3.gb')
    parser.add_argument('--output', help='the RNA file to edit', type=str,
                        default='./coffee1_1.jpg')
    args = parser.parse_args()
    mydrawgenemap(args.input, args.output)
