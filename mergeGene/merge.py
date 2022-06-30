import os

def mergeFile():
    # 以只读模式打开inputOriginally.fasta文件，并读取文件内容
    with open('inputOriginally.fasta', 'r') as OriginallyMerge:
        data = OriginallyMerge.readlines()

        # 原基因物种名称
        name = data[0][1:]
        # 原基因序列
        gene = data[1]

    with open('inputContrast.fasta', 'r') as ContrasMerge:
        # 按行读取文件内容
        data = ContrasMerge.readlines()
        # 计数器
        i = 0
        # 待拼接
        txt = ''
        # 对读取的内容进行拼接
        for line in data:
            if i % 2 == 0 :
                names = name.replace('\n', '') + '-' + data[i][1:]
                txt = txt + names
            else:
                genes = gene + '\n' + data[i]
                txt = txt + genes + '\n'
            i = i + 1
    # 将拼接的文件写入到待计算的axt文件中
    with open("result.axt", "w", encoding='utf-8') as f:
        f.write(str(txt))
        f.close()
    # print(txt)

# 程序入口
if __name__ == '__main__':
    mergeFile()
    os.system('.\KaKs.exe -i result.axt -o txt.kaks -c 5 -m NG')
    kaks =[]
    with open('txt.kaks', 'r') as f:
        # 按行读取文件内容
        data = f.read().split()
        print(data)
        for i in range(4, len(data)):
            if (i - 4) % 22 == 0:
                kaks.append(data[i])
    print(kaks)