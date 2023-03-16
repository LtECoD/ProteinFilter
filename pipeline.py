import os
import sys
import time
import math
import gzip
import shutil
import tarfile
import argparse
from Bio.PDB import FastMMCIFParser
from multiprocessing import Pool


def check(structure, args):
    chain = structure[0]['A']
    if len(chain) > args.max_len:       # too long
        return 1
    if len(chain) < args.min_len:       # too short
        return 2

    chain.atom_to_internal_coordinates()
    residues = list(chain.get_residues())

    # filter plddt
    for res in residues:
        if res["CA"].bfactor < args.plddt:
            return 3

    # amino acid check
    for idx, res in enumerate(residues):
        if len(res.child_dict) <= 0:    # no atom
            return 4
        if not res.internal_coord.is20AA:   # non-standard aa
            return 4
        if idx < (len(residues)-1):
            peptide_bond = res.internal_coord.get_length("C:1N")
            if (peptide_bond is None) or \
                (peptide_bond > args.max_peptide_bond) or \
                    (peptide_bond < args.min_peptide_bond):
                return 4
    return 0


def filter_func(pid, gzipfps, args):
    parser = FastMMCIFParser()

    valid = 0
    too_long = 0
    too_short = 0
    lower_plddt = 0
    bad_res = 0

    for fp in gzipfps:
        gzipf = gzip.open(fp, "rt")
        name = os.path.basename(fp)[:-7]
        structure = parser.get_structure(structure_id=name, filename=gzipf)

        code = check(structure, args)
        if code == 0:
            valid += 1
            shutil.move(fp, os.path.join(args.result_dir, os.path.basename(fp)))
        elif code == 1:
            too_long += 1
        elif code == 2:
            too_short += 1
        elif code == 3:
            lower_plddt += 1
        elif code == 4:
            bad_res += 1
        else:
            raise ValueError(f"Unvalid check code: {code}")

    return [valid, too_long, too_short, lower_plddt, bad_res]
        

def unzip_fun(pid, tar_fp, tmp_dir):
    tarf = tarfile.open(tar_fp)
    cif_mems = [subf for subf in tarf.getmembers() if subf.name.endswith(".cif.gz")]
    tarf.extractall(path=tmp_dir, members=cif_mems)
    return len(cif_mems)


def error_callback(err):
    print(f'Error: {str(err)}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--af_dir", type=str)
    parser.add_argument("--tmp_dir", type=str, default="/lustre/Data/tmp_cif")
    parser.add_argument("--result_dir", type=str, default="/lustre/home/yangsen/database/AlphaFold_filtered")
    parser.add_argument("--log", type=str, default="./log.txt")
    parser.add_argument("--processes", type=int) 

    parser.add_argument("--plddt", type=float)
    parser.add_argument("--min_len", type=int)
    parser.add_argument("--max_len", type=int)
    parser.add_argument("--min_peptide_bond", type=float)
    parser.add_argument("--max_peptide_bond", type=float)
    args = parser.parse_args()

    logger_handler = open(args.log, "w")
    logger_handler.write("CIF\tRatio\tValid\tLong\tShort\tPLDDT\tRES\tExtract\tFilter\n")

    species_tars = os.listdir(args.af_dir)
   
    if os.path.exists(args.result_dir):
        shutil.rmtree(args.result_dir)
    os.makedirs(args.result_dir, exist_ok=True)

    tar_idx = -1
    while True:
        start_time = time.time()
        # clear tmp directory
        tmp_cif_num = 0
        shutil.rmtree(args.tmp_dir)
        os.makedirs(args.tmp_dir)

        # unzip 
        while True:
            unzip_pool = Pool(processes=args.processes)
            unzip_results = []
            for pid in range(args.processes):
                tar_idx = tar_idx + 1
                if tar_idx >= len(species_tars):
                    break
                unzip_results.append(unzip_pool.apply_async(
                    unzip_fun,
                    args=(pid, 
                          os.path.join(args.af_dir, species_tars[tar_idx]),
                          args.tmp_dir),
                    error_callback=error_callback))
                
            unzip_pool.close()
            unzip_pool.join()

            tmp_cif_num += sum([r.get() for r in unzip_results])
            if (tmp_cif_num > args.processes * 100) or (tar_idx >= len(species_tars)):
                break
        extract_time = time.time()

        # filter
        filter_pool = Pool(processes=args.processes)
        gzipfps = [os.path.join(args.tmp_dir, fn) for fn in os.listdir(args.tmp_dir)]
        assert len(gzipfps) == tmp_cif_num
        batch = math.ceil(tmp_cif_num / args.processes)

        results = []
        for pid in range(args.processes):
            l, r = pid * batch, min((pid+1) * batch, tmp_cif_num)
            results.append(filter_pool.apply_async(
                filter_func,
                args=(pid,
                      gzipfps[l: r],
                      args),
                error_callback=error_callback))
        filter_pool.close()
        filter_pool.join()
        filter_time = time.time()

        _tot = 0
        tot_valid, tot_long, tot_short, tot_plddt, tot_res = 0, 0, 0, 0, 0
        for idx in range(len(results)):
            _result = results[idx].get()
            _tot += sum(_result)

            tot_valid += _result[0]
            tot_long += _result[1]
            tot_short += _result[2]
            tot_plddt += _result[3]
            tot_res += _result[4]
        assert _tot == tmp_cif_num
        logger_handler.write(
            f"{tmp_cif_num}\t{round(tot_valid/tmp_cif_num*100, 2)}%" + 
            f"\t{tot_valid}\t{tot_long}\t{tot_short}\t{tot_plddt}\t{tot_res}" +
            f"\t{round((extract_time-start_time)/60, 2)}" +
            f"\t{round((filter_time-extract_time)/60, 2)}\n")
        logger_handler.flush()

        if tar_idx >= len(species_tars):
            break

    logger_handler.close()

