import os
import matplotlib
from matplotlib import pyplot as plt

from flask import Flask, request

app = Flask(__name__)


# app.route装饰器映射URL和执行的函数。这个设置将根URL映射到了hello_world函数上
@app.route('/hello')
def hello_world():
    return 'Hello World!'

# 对fasta文件进行处理拆分并生成axt文件函数
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

    # 对原基因进行切割,并保存在listWindowOrigin里
    listWindowOrigin = []

    for i in range(0, geneLen, int(step)):
        listWindowOrigin.append(gene[i:i+int(size)])


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
                listWindowContrast = []

                for j in range(0, geneLen, int(step)):
                    listWindowContrast.append(data[i][j:j + int(size)])

                    # 拼接名字
                    txt = txt + names

                    # 原基因和需要对比基因进行拼接
                    genes = listWindowOrigin[flag] + '\n' + listWindowContrast[flag]
                    flag = flag + 1
                    txt = txt + genes + '\n\n'
            i = i + 1


    # 将拼接的文件写入到待计算的axt文件中
    with open("result.axt", "w", encoding='utf-8') as f:
        f.write(str(txt))
        f.close()
    # print(txt)
    return geneLen, namelist

# 计算单个文件的kaks函数
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

        # 根据用户选择处理文件，将文件拆分为两个文件
        zeroSpecieStr = ''
        otherSpeciesStr = ''

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

    # 字典保存每一条比对数据
    yKaks = {}

    number = 0
    speciesKaks = []

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

    return yKaks

# 计算指定文件的kaks值，并以字典的数据格式返回
@app.route('/getfasfilekaks', methods=['POST'])
def getFasFileKaks():
    fileName = request.values.get('fileName')

    zeroSpecie = request.values.get('zeroSpecie')

    windowSize = request.values.get('windowSize')

    slidingStep = request.values.get('slidingStep')

    kaks = calculate(fileName, zeroSpecie, windowSize, slidingStep)

    return {
        'msg': '返回指定文件的Kaks值',
        'status': True,
        'result': kaks,
    }

# 计算指定path下的所有文件的kaks值，并以字典的数据格式返回
@app.route('/getallfasfilekaks', methods=['POST'])
def getAllFasFileKaks():

    path = request.values.get('path')

    zeroSpecie = request.values.get('zeroSpecie')

    windowSize = request.values.get('windowSize')

    slidingStep = request.values.get('slidingStep')

    print(path)

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

    # ----------------循环遍历计算fastaFile下所有文件---------
    for f in fileNameList:
        fkaks = calculate(f, zeroSpecie, windowSize, slidingStep)
        kaksList[f] = fkaks

    return {
        'msg': '返回所有15个文件的Kaks值',
        'status': True,
        'result': kaksList,
    }

# 程序入口
if __name__ == '__main__':

    app.run(host='0.0.0.0',port=8520)
