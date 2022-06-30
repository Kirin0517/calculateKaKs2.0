import os
import matplotlib
import numpy as np
from  matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator


def mergeFile(size):

    # 物种名数组，用于后面作图用
    namelist = []

    # 以只读模式打开inputOriginally.fasta文件，并读取文件内容
    with open('inputOriginally.fasta', 'r') as OriginallyMerge:
        data = OriginallyMerge.readlines()

        # 原基因物种名称
        name = data[0][1:]

        # 原基因序列
        gene = data[1]
        # print(gene)

    # 滑窗移动次数
    geneLen = len(gene)-int(size)

    # 对原基因进行切割
    listWindowOrigin = []
    for i in range(0, geneLen):
        listWindowOrigin.append(gene[i:i+int(size)])
    # print(listWindowOrigin)
    # print(len(listWindowOrigin))


    with open('inputContrast.fasta', 'r') as ContrasMerge:
        # 按行读取文件内容
        data = ContrasMerge.readlines()
        # 计数器
        i = 0
        # 待拼接axt输入数据
        txt = ''


        # 对读取的内容进行拼接
        for line in data:
            # 物种名称拼接部分
            if i % 2 == 0:
                names = name.replace('\n', '') + '-' + data[i][1:]
                namelist.append(names)
            # 基因拼接部分
            else:
                # print(data[i])
                listWindowContrast = []
                # print(data[i])
                # print(len(data[i]))

                for j in range(0, geneLen):
                    # print(j)
                    listWindowContrast.append(data[i][j:j + int(size)])

                    # 拼接名字
                    txt = txt + names
                    # print(names)
                    # 原基因和需要对比基因进行拼接
                    genes = listWindowOrigin[j] + '\n' + listWindowContrast[j]
                    # print(genes)
                    txt = txt + genes + '\n\n'
                # print(j)
            i = i + 1


    # 将拼接的文件写入到待计算的axt文件中
    with open("result.axt", "w", encoding='utf-8') as f:
        f.write(str(txt))
        f.close()
    # print(txt)
    return geneLen, namelist

def drawKaKs(speciesList, xAxisDate, yKaksDates):
    # print(speciesList, xAxisDate, yKaksDate)
    plt.figure(figsize=(6, 4))  # 创建绘图对象

    # 是否显示背景网格
    matplotlib.rc('axes', grid=False)

    x = []
    for i in range(0, xAxisDate):
         x.append(i)
    print(x)
    print(yKaksDates)
    # x = [0, 1, 2, 3, 4, 5, 6]
    # y = [0.3, 0.4, 2, 5, 3, 4.5, 4]
    yKaks = []
    for i in range(0, len(yKaksDates)):

        # yKaksDates[i]数据是文本，需要格式转换，不然图像y坐标轴不排序！！！！！！！！！！！！！！！
        if yKaksDates[i] == 'NA':
            yKaks.append(0)
        else:
            yKaks.append(float(yKaksDates[i]))
        if (i + 1) % len(x) == 0:
            print(yKaks)
            plt.plot(x, yKaks , linewidth=0.3)  # 在当前绘图对象绘图（X轴，Y轴，蓝色虚线，线宽度）
            yKaks =[]

    plt.xlabel("sliding window (staring codon)")  # X轴标签
    plt.ylabel("Ka / Ks")  # Y轴标签
    plt.title("Line plot")  # 图标题


    plt.show()  # 显示图
    plt.savefig("line.jpg")  # 保存

# 程序入口
if __name__ == '__main__':

    # 设置滑窗大小
    windowSize = input('设置滑窗大小：')
    param = mergeFile(windowSize)

    xAxisLen = param[0]
    speciesNameList = param[1]
    # print(data)
    os.system('.\KaKs.exe -i result.axt -o txt.kaks -c 5 -m NG')
    kaks =[]
    with open('txt.kaks', 'r') as f:
        # 按行读取文件内容
        data = f.read().split()
        # print(data)
        for i in range(26, len(data)):
            if (i - 4) % 22 == 0:
                kaks.append(data[i])
    print(kaks)
    # drawKaKs(speciesNameList, xAxisLen, kaks)
