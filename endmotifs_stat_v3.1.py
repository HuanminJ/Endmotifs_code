import pandas as pd
import os, sys,argparse
import subprocess
import glob
import pysam
from reverse_complementary_DNA import complementary
from functools import reduce
from motifs_generate import motifs_GC
from scipy.stats import entropy
from ref_base_stat1 import run_base_stat
from myargs_parse import get_args



def stat_endmotifs(indir,
                   sample_list,
                   ref,
                   num_base,
                   suffix,
                   outdir,
                   std_AT=None,
                   std_GC=None,
                   start=None,
                   end=None,
                   vref=False):
    n = num_base
    motif_set_n = motifs_GC(n=num_base)
    f = open(os.path.join(indir, sample_list), "r")
    samples = []
    list_total = []
    for s in f.readlines():
        s = s.strip()
        print(s)
        item = os.path.join(indir, s, s + suffix)
        seqDic_t = dict((k, []) for k in range(1, n + 1))
        if os.path.getsize(item) > 1:
            samples.append(s)
            sam_ali = pysam.AlignmentFile(item, 'rb')  #bam文件打开
            for read in sam_ali.fetch(reference=ref, start=start, stop=end):
                qaseq = read.query_alignment_sequence
                flag = read.flag
                ref = read.reference_name
                if vref == False:
                    if ref == None:
                        if flag == 83:
                            for k in range(1, n + 1):
                                seqDic_t[k].append(complementary(qaseq[-k:]))
                        elif flag == 99:
                            for k in range(1, n + 1):
                                seqDic_t[k].append(qaseq[:k])
                        else:
                            pass
                    else:
                        if flag == 83 and ref == ref:
                            for k in range(1, n + 1):
                                seqDic_t[k].append(complementary(qaseq[-k:]))
                        elif flag == 99 and ref == ref:
                            for k in range(1, n + 1):
                                seqDic_t[k].append(qaseq[:k])
                        else:
                            pass
                elif vref == True and ref != None:
                    if flag == 83 and ref != ref:
                        for k in range(1, n + 1):
                            seqDic_t[k].append(complementary(qaseq[-k:]))
                    elif flag == 99 and ref != ref:
                        for k in range(1, n + 1):
                            seqDic_t[k].append(qaseq[:k])
                    else:
                        pass
                else:
                    print(
                        "Parameter -v/--vref must be combined with -r/--ref,please check again!"
                    )
            sam_ali.close()

            #统计末端相应末端碱基类型reads数量
            total_count = {
                m: seqDic_t[k].count(m)
                for k in range(1, n + 1) for m in motif_set_n[k]
            }
            end_outpath = os.path.join(indir, outdir, "samples_data")
            os.makedirs(end_outpath, exist_ok=True)
            df_t = pd.DataFrame(total_count, index=[0])
            list_total.append(df_t)
            df_t.to_csv(
                os.path.join(end_outpath, f"{s}_endmotifs_count.csv"),
                index=False)
        else:
            print(f"{s} bam too small and skipped!")

    #末端结果合并
    df_total = pd.concat(list_total, axis=0)
    df_total = df_total.reset_index()
    df_total = df_total.drop('index', axis=1)
    dfs = pd.concat([pd.DataFrame({"sample": samples}), df_total], axis=1)
    dfs.to_csv(
        os.path.join(indir, outdir, f"total_endmotifs_count.csv"), index=False)
    f.close()

    for i in range(1, num_base + 1):
        df_total.loc[:,
                     motif_set_n[i]] = df_total.loc[:, motif_set_n[i]].apply(
                         lambda x: x / x.sum(), axis=1)
        #compute endmotifs diversity
        pk = df_total.loc[:, motif_set_n[i]].T
        centropy = entropy(pk, base=len(motif_set_n[i]))
        mds = pd.DataFrame({"sample": samples, "centropy": centropy})
        mds.to_csv(
            os.path.join(outdir, f"endmotifs_diversity_{i}.csv"), index=False)
    if std_AT and std_GC:
        A = pd.DataFrame({"A_std": df_total['A'] / float(std_AT)})
        T = pd.DataFrame({"T_std": df_total['T'] / float(std_AT)})
        G = pd.DataFrame({"G_std": df_total['G'] / float(std_GC)})
        C = pd.DataFrame({"C_std": df_total['C'] / float(std_GC)})
        std_atgc = pd.concat([pd.DataFrame({"sample": samples}), A, T, G, C],axis=1)
        std_atgc.to_csv(
        os.path.join(outdir, "one_base_endmotifs_preference.csv"), index=False)
    df_total = pd.concat([pd.DataFrame({"sample": samples}), df_total], axis=1)
    df_total.to_csv(
        os.path.join(outdir, "total_endmotifs_ratio.csv"), index=False)
    print('processed done,file saved!')

    return None


if __name__ == "__main__":
    args = get_args()
    region = args.start_end
    if args.fasta != None and args.ref != None:
        std_AT, std_GC = run_base_stat(
            fasta=args.fasta,
            outdir=args.outdir,
            ref=args.ref,
            start=region[0],
            end=region[1],
            vref=args.vref)
        stat_endmotifs(
            indir=args.indir,
            sample_list=args.samples,
            ref=args.ref,
            num_base=args.num_base,
            suffix=args.suffix,
            outdir=args.outdir,
            std_AT=std_AT,
            std_GC=std_GC,
            start=region[0],
            end=region[1],
            vref=args.vref)
    else:
        stat_endmotifs(
            indir=args.indir,
            sample_list=args.samples,
            ref=args.ref,
            num_base=args.num_base,
            suffix=args.suffix,
            outdir=args.outdir,
            start=region[0],
            end=region[1],
            vref=args.vref)