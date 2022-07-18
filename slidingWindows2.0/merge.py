import os
import matplotlib
from matplotlib import pyplot as plt

from flask import Flask, request

app = Flask(__name__)


# app.route装饰器映射URL和执行的函数。这个设置将根URL映射到了hello_world函数上
@app.route('/hello')
def hello_world():
    return 'Hello World!'

# 文件拼接函数
def mergeFile(size, step):

    # 物种名数组，用于后面作图用
    namelist = []

    # 以只读模式打开inputOriginally.fasta文件，并读取文件内容
    with open('./compareFasta/zeroSpecieFile.fasta', 'r') as OriginallyMerge:
        data = OriginallyMerge.readlines()

        # 原基因物种名称
        name = data[0][1:]

        # 原基因序列
        gene = data[1]

    # 滑窗需要移动的长度
    geneLen = len(gene)-int(size)+int(step)
    # print(geneLen)
    # 对原基因进行切割,并保存在listWindowOrigin里
    listWindowOrigin = []

    for i in range(0, geneLen, int(step)):
        listWindowOrigin.append(gene[i:i+int(size)])
        # print(i)

    # print(listWindowOrigin)
    # print(len(listWindowOrigin))


    with open('./compareFasta/otherSpeciesFile.fasta', 'r') as ContrasMerge:

        # 按行读取文件内容
        data = ContrasMerge.readlines()

        # 计数器
        i = 0

        # 待拼接axt输入数据
        txt = ''

        # 对读取的内容进行拼接
        for line in data:

            # 计数器
            flag = 0

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

                for j in range(0, geneLen, int(step)):
                    # print(j)
                    listWindowContrast.append(data[i][j:j + int(size)])
                    # print(listWindowContrast)
                    # 拼接名字
                    txt = txt + names
                    # print(names)
                    # 原基因和需要对比基因进行拼接
                    genes = listWindowOrigin[flag] + '\n' + listWindowContrast[flag]
                    flag = flag + 1
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


# 画图函数
def drawKaKs(speciesList, xAxisDate, yKaksDates, slidingStep):
    # print(speciesList, xAxisDate, yKaksDates)

    # print(len(speciesList), speciesList)
    print(yKaksDates)
    plt.figure(figsize=(6, 4))  # 创建绘图对象

    # 是否显示背景网格
    matplotlib.rc('axes', grid=False)

    # x轴数据赋值
    x = []
    number = 1

    for i in range(0, xAxisDate):
         x.append(i)
    print(x)
    # print(yKaksDates)

    yKaks = []
    for i in range(0, len(yKaksDates)):

        # yKaksDates[i]数据是文本，需要格式转换，不然图像y坐标轴不排序！！！！！！！！！！！！！！！
        if yKaksDates[i] == 'NA':
            yKaks.append(0)
        else:
            yKaks.append(float(yKaksDates[i]))
        if (i + 1) % (len(x) / int(slidingStep)) == 0:
            print('第'+str(number) + '次画图： ', yKaks)
            plt.plot(x, yKaks, linewidth=0.3)  # 在当前绘图对象绘图（X轴，Y轴，蓝色虚线，线宽度）
            plt.legend(speciesList[number])
            print(speciesList[number])
            number = number + 1
            yKaks =[]

    plt.xlabel("sliding window (staring codon)")  # X轴标签
    plt.ylabel("Ka / Ks")  # Y轴标签
    plt.title("Line plot")  # 图标题

    plt.savefig("line.jpg")  # 保存
    plt.show()  # 显示图


def calculate(inputfile, zeroSpecie, windowSize, slidingStep):
    # 打开文件
    with open(inputfile, 'r') as f:

        # 按行读取文件
        data = f.readlines()

        # 保存文件下所有物种名称
        speciesNameList = []

        # 计数器，用于筛选出文件中的所有物种名
        k = 0
        for line in data:  # 遍历文件，一行行遍历，读取文本
            if k % 2 == 0:
                speciesNameList.append(line[1:])
            k += 1

        # print(speciesNameList)
        # print('Please select species zero from above species:')

        # 假设用户选择第三个物种作为0物种，其余为将要比对的物种
        # zeroSpecie = speciesNameList[2]

        # 根据用户选择处理文件，将文件拆分为两个文件
        zeroSpecieStr = ''
        otherSpeciesStr = ''

        # print(data)
        # 遍历文件，一行行遍历，读取文本
        # flag用于标志读取两行
        flag = 0
        for line in data:
            if line[1:].strip() == zeroSpecie:          # 由于line包含换行，所以.strip()去除换行符再进行比较
                zeroSpecieStr = zeroSpecieStr + line
                flag = 1
            elif flag == 1:
                zeroSpecieStr = zeroSpecieStr + line.strip()    # 大坑！！！！！！！！！！！！！！！！！
                flag = 0
            else:
                otherSpeciesStr = otherSpeciesStr + line

        # 将拼接的文件写入到待计算的fasta文件中
        with open("./compareFasta/zeroSpecieFile.fasta", "w", encoding='utf-8') as f:
            print(zeroSpecieStr)
            f.write(str(zeroSpecieStr))
            f.close()

        with open("./compareFasta/otherSpeciesFile.fasta", "w", encoding='utf-8') as f:
            print(otherSpeciesStr)
            f.write(str(otherSpeciesStr))
            f.close()

    # ---------数据计算-------------------
    param = mergeFile(windowSize, slidingStep)
    # param = mergeFile(12, 3)
    # 每个物种对比次数
    xAxisLen = param[0]
    print(xAxisLen)
    # 保存物种名称
    speciesNameList = param[1]
    os.system('.\KaKs.exe -i result.axt -o txt.kaks -c 5 -m NG')

    print(speciesNameList)

    kaks = []
    with open('txt.kaks', 'r') as f:
        # 按行读取文件内容
        data = f.read().split()
        # print(data)
        for i in range(26, len(data)):
            if (i - 4) % 22 == 0:
                kaks.append(data[i])
    print(kaks)

    # 字典保存每一条对比数据
    yKaks = {}

    number = 0
    speciesKaks = []
    # print(xAxisLen)
    # print(len(kaks))
    for i in range(0, len(kaks)):

        # yKaksDates[i]数据是文本，需要格式转换，不然图像y坐标轴不排序！！！！！！！！！！！！！！！
        if kaks[i] == 'NA' or kaks[i] == '-nan(ind)':
            speciesKaks.append(0)
        else:
            speciesKaks.append(float(kaks[i]))
        if (i + 1) % (xAxisLen / int(slidingStep)) == 0:
            yKaks[speciesNameList[number]] = speciesKaks

            print(speciesNameList[number])
            print(speciesKaks)
            number = number + 1
            speciesKaks = []

    # print(len(yKaks))
    # print(yKaks)
    return yKaks


def calculateAllFasFile(path, zeroSpecie, windowSize, slidingStep):

    # 读取fastafile文件夹下所有文件
    fastafile = os.listdir(path)

    # 保存文件夹名到filelist变量中
    fileNameList = []

    # 保存各个文件计算的kaks
    kaksList = {}

    for filename in fastafile:  # 遍历文件夹
        if not os.path.isdir(filename):  # 判断是否是文件夹，不是文件夹才打开
            if filename.endswith('.fas'):
                fileNameList.append(path + filename)

    print(fileNameList)

    # ----------------循环遍历计算fastaFile下所有文件---------
    for f in fileNameList:
        fkaks = calculate(f, zeroSpecie, windowSize, slidingStep)
        kaksList[f] = fkaks

    return kaksList

# 程序入口
if __name__ == '__main__':

    # app.run(host='0.0.0.0',port=8520)

    # path = './fastaFile/'
    #
    # # 读取fastafile文件夹下所有文件
    # fastafile = os.listdir(path)
    #
    # # 保存文件夹名到filelist变量中
    # fileNameList = []
    #
    # for filename in fastafile:  # 遍历文件夹
    #     if not os.path.isdir(filename):  # 判断是否是文件夹，不是文件夹才打开
    #         if filename.endswith('.fas'):
    #             fileNameList.append(path + filename)
    #
    # print(fileNameList)
    # #
    # #
    # # calculate('./fastaFile/cox1_mafft.fas', '10_Nesodiprion_japonicus_-_Atp6_Test_sp._mitochondrion', 66, 3)
    #
    # # ----------------循环遍历计算fastaFile下所有文件---------
    # for f in fileNameList:
    #     fkaks = calculate(f, 'Dimorphopteryx_sp__1_GYN_2021_OK329944', 66, 3)
    #     print(fkaks)

    # ----------------测试---------------------
    # 设置滑窗大小
    # windowSize = input('Enter the slider size (requires a multiple of the window size to 3)：')
    # while int(windowSize) % 3 != 0:
    #     print('Enter the data is not a multiple of 3, please re-enter!')
    #     windowSize = input('Enter the slider size (requires a multiple of the window size to 3)：')
    #
    # slidingStep = input('Sets the sliding window step size：')
    # # param = mergeFile(windowSize)
    #
    # param = mergeFile(windowSize, slidingStep)
    # # param = mergeFile(12, 3)
    # # 每个物种对比次数
    # xAxisLen = param[0]
    # # print(xAxisLen)
    # # 保存物种名称
    # speciesNameList = param[1]
    # os.system('.\KaKs.exe -i result.axt -o txt.kaks -c 5 -m NG')
    # kaks =[]
    # with open('txt.kaks', 'r') as f:
    #     # 按行读取文件内容
    #     data = f.read().split()
    #     # print(data)
    #     for i in range(26, len(data)):
    #         if (i - 4) % 22 == 0:
    #             kaks.append(data[i])
    # print(kaks)
    #
    # # 字典保存每一条对比数据
    # yKaks = {}
    #
    # number = 0
    # speciesKaks = []
    # print(xAxisLen)
    # print(len(kaks))
    # for i in range(0, len(kaks)):
    #
    #     # yKaksDates[i]数据是文本，需要格式转换，不然图像y坐标轴不排序！！！！！！！！！！！！！！！
    #     if kaks[i] == 'NA':
    #         speciesKaks.append(0)
    #     else:
    #         speciesKaks.append(float(kaks[i]))
    #     if (i + 1) % (xAxisLen / int(slidingStep)) == 0:
    #         yKaks[speciesNameList[number]] = speciesKaks
    #         print(speciesNameList[number])
    #         print(speciesKaks)
    #         number = number + 1
    #         speciesKaks = []
    #
    # print(len(yKaks))
    # print(yKaks)
    # drawKaKs(speciesNameList, xAxisLen, kaks, slidingStep)

    x = calculateAllFasFile('./fastaFile/', 'Dimorphopteryx_sp__1_GYN_2021_OK329944', 66, 3)
    print(x)
    print(len(x))
    print(x['./fastaFile/rrnS_mafft.fas'])