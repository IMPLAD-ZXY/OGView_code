import os
from PIL import Image, ImageDraw, ImageFont
import numpy as np
from PIL import Image
from PIL import ImageDraw, ImageFont

data = {}
width = 3000 # 原图宽度
width_ratio = width / 10000
height = int(5000 * width_ratio)  # 原图高度
left_padding = int(1450 * width_ratio)  # 左侧填充大小
right_padding = int(400 * width_ratio)
fontpath1 = "/home/xyzhang/pmgdraw/bin/arial.ttf"
fontpath2 = "/home/xyzhang/pmgdraw/bin/arialbd.ttf"
fontpath3 = "/home/xyzhang/pmgdraw/bin/arialbi.ttf"
font_mRNA = int(220 * width_ratio)
tail_height = int(300 * width_ratio)  # 最底部的白边


class gene_img:
    def __init__(self):
        self.dna2 = []
        self.dna3 = None
        self.name = None
        return

    def add_dna2(self, img):
        self.dna2.append(img)

    def set_dna3(self, img):
        self.dna3 = img

    def set_name(self, img):
        self.name = img

    def splice(self):
        img_list = [self.name] + self.dna2 + [self.dna3]
        img = np.concatenate(img_list, axis=0)
        shape = img.shape
        img = Image.fromarray(img)
        ImageDraw.Draw(img).rounded_rectangle(
            ((int(100 * width_ratio), int(700 * width_ratio)), (shape[1] - 100, shape[0] - 0)),
            radius=int(150 * width_ratio),
            width=int(30 * width_ratio), outline=(28, 187, 140))  # 给每个基因外边加的框 右 上 左 下
        return img


# 从dna3中提取所有基因名字
def obtain_gene_names(dir_name: str, flag, title):
    """
    提取input2中所有的基因名字
    加载input2中所有的图片
    添加表头tran，cis
    """
    data.clear()
    file_names = os.listdir(dir_name)
    for file_name in file_names:
        gene_name = file_name[:-4]  # 文件命名方式为'基因名.jpg'
        data[gene_name] = gene_img()
    header_img = Image.open(os.path.join(dir_name, file_names[0])).resize((width, height), Image.ANTIALIAS)

    header_img = np.pad(header_img.crop((int(30 * width_ratio), 0, int(9998 * width_ratio), int(850 * width_ratio))), (
    (0, 0),
    (int((left_padding + right_padding) / 2), left_padding + right_padding - int((left_padding + right_padding) / 2)),
    (0, 0)), 'constant',
                        constant_values=255)
    flag_img = Image.fromarray(
        np.ones((int(1400 * width_ratio), header_img.shape[1], 3), dtype=np.uint8) * 255)  # 最上方添加空白距离
    # 原先坐标是文本框左上角坐标，利用‘anchor='mm'’将坐标改为文本框中心坐标，实现将名字放在中间
    font1 = ImageFont.truetype(fontpath1, int(250 * width_ratio))
    font2 = ImageFont.truetype(fontpath3, int(350 * width_ratio))
    ImageDraw.Draw(flag_img).text((int(header_img.shape[1] / 2), int(1200 * width_ratio)), flag + '-splicing genes',
                                  font=font1, fill=(0, 0, 0), anchor='mm')
    ImageDraw.Draw(flag_img).text((int(header_img.shape[1] / 2), int(900 * width_ratio)), 'Plastid Genome',
                                  font=font1, fill=(0, 0, 0), anchor='mm')
    ImageDraw.Draw(flag_img).text((int(header_img.shape[1] / 2), int(400 * width_ratio)), title, font=font2,
                                  fill=(0, 0, 0), anchor='mm')  # 表头物种名位置
    # header_img = np.concatenate([header_img, flag_img])
    return flag_img


def load_input1_fig(dir_name: str, CHR: str):
    """
    加载input1中所有的图片

    """
    for gene_name in data.keys():
        if os.path.exists(os.path.join(dir_name, gene_name + '+.jpg')) and os.path.exists(
                os.path.join(dir_name, gene_name + '-.jpg')):
            p_img = Image.open(os.path.join(dir_name, gene_name + '+.jpg')).resize((width, height), Image.ANTIALIAS)
            n_img = Image.open(os.path.join(dir_name, gene_name + '-.jpg')).resize((width, height), Image.ANTIALIAS)
            data[gene_name].add_dna2(splice_p_n_imgs(p_img, n_img, CHR))


def splice_p_n_imgs(p_img, n_img, CHR):
    """
    裁剪input1的图片，并将同一染色体的同一基因正负链拼接
    给裁剪后的图片左右两边扩充空白
    在左右两边加5‘+和3’-的文字
    """
    # 裁剪n_img
    n_img = n_img.crop(
        (int(30 * width_ratio), int(2420 * width_ratio), int(9998 * width_ratio), int(3480 * width_ratio)))

    # 裁剪p_img
    p_img = p_img.crop(
        (int(30 * width_ratio), int(1510 * width_ratio), int(9998 * width_ratio), int(2570 * width_ratio)))

    # 拼接
    img = np.concatenate((p_img, n_img), axis=0)  # 纵向拼接

    # 边缘扩充，左边扩充。。。像素，右边扩充。。。像素
    img1 = np.pad(img, ((0, int(20 * width_ratio)), (left_padding, right_padding), (0, 0)), 'constant',
                  constant_values=255)
    img1 = Image.fromarray(img1)

    # 加字
    font = ImageFont.truetype(fontpath2, font_mRNA)
    right_index = left_padding + int(9800 * width_ratio)
    left_index = left_padding
    high_index = int(1070 * width_ratio)
    low_index = int(850 * width_ratio)
    ImageDraw.Draw(img1).text((left_index, high_index), "3'", font=font, fill=(0, 0, 0))
    ImageDraw.Draw(img1).text((left_index, low_index), "5'", font=font, fill=(0, 0, 0))
    ImageDraw.Draw(img1).text((right_index, low_index), "3'", font=font, fill=(0, 0, 0))
    ImageDraw.Draw(img1).text((right_index, high_index), "5'", font=font, fill=(0, 0, 0))
    CHR = '\n'.join([CHR[i:i+9] for i in range(0, len(CHR), 9)])
    ImageDraw.Draw(img1).text((left_padding - int(1250 * width_ratio), int((high_index + low_index) / 2)), CHR,
                              font=font, fill=(0, 0, 0))

    return img1


def load_input2_fig(dir_name: str):
    """
    裁剪input2的图片
    给裁剪后的图片左右两边扩充空白
    在左边加mRNA的文字
    添加每个基因的名字，加粗倾斜
    """
    for gene_name in data.keys():
        img = Image.open(os.path.join(dir_name, gene_name + '.jpg')).resize((width, height), Image.ANTIALIAS)
        img1 = np.pad(img.crop(
            (int(30 * width_ratio), int(1700 * width_ratio), int(9998 * width_ratio), int(2650 * width_ratio))),
                      ((0, 0), (left_padding, right_padding), (0, 0)), 'constant',
                      constant_values=255)
        img_shape = img1.shape
        img1 = Image.fromarray(img1)

        font = ImageFont.truetype(fontpath2, font_mRNA)
        ImageDraw.Draw(img1).text(
            (left_padding - int(1000 * width_ratio) + int(100 * width_ratio), int(550 * width_ratio)), "mRNA",
            font=font, fill=(0, 0, 0))
        data[gene_name].set_dna3(img1)

        name_img = Image.fromarray(np.ones((int(800 * width_ratio), img_shape[1], 3), dtype=np.uint8) * 255)
        # 原先坐标是文本框左上角坐标，利用‘anchor='mm'’将坐标改为文本框中心坐标，实现将名字放在中间
        font = ImageFont.truetype(fontpath3, int(300 * width_ratio))  # 每个基因名字的字体大小
        ImageDraw.Draw(name_img).text((int(img_shape[1] / 2), int(500 * width_ratio)), gene_name.replace(' ', ''), font=font,
                                      fill=(0, 0, 0), anchor='mm')
        data[gene_name].set_name(name_img)


def splice_genes(header_img):
    """
    给整个图片添加表头
    表头包括两行，物种名和mitochondrial genome
    把所有基因和图头拼在一起
    """
    img_list = [header_img]
    for gene_name in data.keys():
        img_list.append(data[gene_name].splice())
    img_list.append(np.ones((tail_height, img_list[0].width, 3), dtype=np.uint8) * 255)
    return np.concatenate(img_list, axis=0)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--input1', help='the RNA file to edit', type=str, nargs='+',
                        default=['D:/Users/Administrator/PycharmProjects/pythonProject/新建文件夹/dna2_1/trans',
                                 'D:/Users/Administrator/PycharmProjects/pythonProject/新建文件夹/dna2_2/trans'])
    parser.add_argument('--input2', help='the RNA file to edit', type=str,
                        default='D:/Users/Administrator/PycharmProjects/pythonProject/新建文件夹/dna3/trans')
    parser.add_argument('--title', help='title', type=str,
                        default='mitochondrion Saccharum officinarum')
    parser.add_argument('--CHR', help='ACCESSION number', type=str, nargs='+',
                        default=['NC_016005', 'NC_016005'])
    parser.add_argument('--output', help='the path to save result', type=str,
                        default='D:/Users/Administrator/PycharmProjects/pythonProject/新建文件夹/trans.jpg')
    parser.add_argument('--flag', help='the path to save result', type=str, default='trans')
    args = parser.parse_args()
    header_img = obtain_gene_names(args.input2, args.flag, args.title)
    load_input2_fig(args.input2)
    for i in range(len(args.input1)):
        load_input1_fig(args.input1[i], CHR=f'chr{i + 1}')
    img = Image.fromarray(splice_genes(header_img))
    img.save(args.output)
