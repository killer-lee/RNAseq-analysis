#####Import Module#####
import logging #记录日志
import sys #提供对一些pyhon解释器密切相关的函数的访问
import os #操作系统交互
import math #数学运算
import time #提供时间函数
import argparse #用于命令行参数解析
import glob #查找符合规则的文件名

#####Description####
usage = '''
@Date    : 2024-11-22 20:30:16
@Author  : kang (kangjichang@outlook.com)
@Link    : None
@Version : None
@FileName: RNA_seq.py
Description:
    A script for upstream analysis of rnaseq. If you want to use it, you need to replace the software installation path I specified. 
    At the same time, it has to be mentioned that for my script, my transcriptome assembly step did not actually work. You can delete this part.
Example:
    python {} [-i input] [-o output]
Step:

'''.format(__file__[__file__.rfind(os.sep) + 1:])


#####HelpFormat#####
class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass 
#argparse.RawDescriptionHelpFormatter允许在帮助信息中保留换行和格式。
#argparse.ArgumentDefaultsHelpFormatter在帮助信息中自动显示参数的默认值。

def show_info(text):
    now_time = time.time()
    logging.info(text)
    return now_time

def run_cmd(cmd):
    logging.info(cmd)
    flag = os.system(cmd)
    if flag != 0:
        logging.error("Command fail: " + cmd)
        exit(2)
    return 0

def run_time(start_time):
    spend_time = time.time() - start_time
    logging.info("Total  spend time : " + fmt_time(spend_time))
    return 0

def fmt_time(spend_time):
    spend_time=int(spend_time)
    day = 24 * 60 * 60
    hour = 60 * 60
    min = 60
    if spend_time < 60:
        return "%ds" % math.ceil(spend_time)
    elif spend_time > day:
        days = divmod(spend_time,day)
        return "%dd%s" % (int(days[0]),fmt_time(days[1]))
    elif spend_time > hour:
        hours = divmod(spend_time,hour)
        return "%dh%s" % (int(hours[0]),fmt_time(hours[1]))
    else:
        mins = divmod(spend_time,min)
        return "%dm%s" % (int(mins[0]),fmt_time(mins[1]))

#####main#########
def main():
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter, description=usage)
    parser.add_argument(
        "-t","--thread",help="thread",dest="thread",type=str)
    parser.add_argument(
        "-r1","--r1",help="r1",dest="r1",type=str)
    parser.add_argument(
        "-r2","--r2",help="r2",dest="r2",type=str)
    parser.add_argument(
        "-r","--ref",help="ref",dest="ref",type=str)
    parser.add_argument(
        "-a","--anno",help="annotation",dest="anno",type=str)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s [line:%(lineno)d][%(levelname)s:] %(message)s',
                        datefmt='%Y-%m-%d  %H:%M:%S',
                        filename="." + "/chipseq.log",
                        filemode='w')
    #####fastqc#########
    folder_name = "FASTQC"
    os.makedirs(folder_name,exist_ok=True)
    fastqc = "/data1/kang/software/FastQC/fastqc"
    run_start = show_info("=======Filter step1:fastqc is start=======")
    cmd = "%s -f fastq -o FASTQC %s %s" % (fastqc, args.r1, args.r2)
    run_cmd(cmd)
    FASTQC=os.path.abspath("FASTQC")
    run_time(run_start)

    #####merge reports#########
    file_name = "FASTQC_Resport"
    os.makedirs(file_name,exist_ok=True)
    multiqc = "/data1/kang/bin/multiqc"
    run_start = show_info("=======Filter step2:merge report is start=======")
    cmd = "%s %s -o FASTQC_Resport" % (multiqc, FASTQC) 
    run_cmd(cmd)
    run_time(run_start)

    #####delet adapter#########
    Trimmomatic = "/data1/kang/software/Trimmomatic-0.39/trimmomatic-0.39.jar"
    name1 = os.path.basename(args.r1).split(".")[0]
    name2 = os.path.basename(args.r2).split(".")[0]
    run_start = show_info("=======Filter step3:delet adapter is start=======")
    cmd = "java -jar %s PE -phred33 -threads %s %s %s %s_clean_data.fastq.gz %s__unpaired_clean.fastq.gz %s__clean_data.fastq.gz %s_unpaired_clean.fastq.gz ILLUMINACLIP:%s/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 HEADCROP:10 SLIDINGWINDOW:4:15 MINLEN:36" % (
    Trimmomatic, args.thread, args.r1, args.r2, name1, name1, name2, name2, Trimmomatic
    )
    run_cmd(cmd)
    run_time(run_start)

    #####mapping#########
    hisat2 = "/data1/kang/software/hisat2-2.2.1/hisat2"
    basename = os.path.basename(args.r1).split("_")[0]
    run_start = show_info("=======Filter step4:mapping is start=======")
    cmd = "%s -t --dta -p %s -x %s -1 %s_clean_data.fastq.gz -2 %s__clean_data.fastq.gz -S %s.sam" % (hisat2, args.thread, args.ref, name1, name2, basename)
    run_cmd(cmd)
    run_time(run_start)

    #####change for########

    run_start = show_info("=======Filter step5:change for is start=======")
    cmd = "samtools view -bS %s.sam > %s.bam" % (basename, basename)
    run_cmd(cmd)
    run_time(run_start)

    #####sort#######
    run_start = show_info("=======Filter step5:Sort is start=======")
    cmd = "samtools sort -o %s.sort.bam %s.bam" % (basename, basename)
    run_cmd(cmd)
    run_time(run_start)

    #####samtools index#######
    run_start = show_info("=======Filter step6:built index is start=======")
    cmd = "samtools index %s.sort.bam" % (basename)
    run_cmd(cmd)
    run_time(run_start)

    #####transcriptome assembly#######
    stringtie = "/data1/kang/stringtie-2.2.3.Linux_x86_64/stringtie"
    run_start = show_info("=======Filter step7:transcriptome asseembly is start=======")
    cmd = "%s -p %s -G %s -o %s.gtf %s.sort.bam" % (stringtie, args.thread, args.anno, basename, basename)
    run_cmd(cmd)
    run_time(run_start)

    #####HTseq counting#######
    htseq = "/data1/kang/anaconda3/bin/htseq-count"
    run_start = show_info("=======Filter step8:HTseq counting is start=======")
    cmd = "%s -f bam %s.sort.bam %s > %s_count.txt" % (htseq, basename, args.anno, basename)
    run_cmd(cmd)
    run_time(run_start)

if __name__ == "__main__":
    try:
        t1 = time.time()
        time1 = time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(t1))
        print("start at:" + time1)

        main()
        t2 = time.time()
        time2 = time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(t2))
        print("End at:" + time2)
        t3 = t2 - t1
        print("Spend time" + str(t3))
    except KeyboardInterrupt:
        sys.stderr.write("Bye Honey\n")
        sys.exit(0)
