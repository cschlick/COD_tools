import argparse
import tqdm
from pathlib import Path
import multiprocessing
from multiprocessing import Pool

from cod_tools.CODFile import CODFile

parser = argparse.ArgumentParser(
                    prog = 'COD filter',
                    description = 'Filter a folder of COD files',
                    epilog = '')
parser.add_argument('--folder', type=str,help="Path to COD folder")


args = parser.parse_args()


cod_dir = Path(args.folder)
cif_files = [path for path in cod_dir.glob("**/*") if path.suffix == ".cif"]
print("Number cif files:",len(cif_files))

work = cif_files#[:10000]
def worker(file):
  return CODFile(file,cod_id=file.stem)

nproc=64
disable_progress = False
results = []
with Pool(nproc) as pool:
  for result in tqdm.tqdm(pool.imap(worker, work),
                          total=len(work),
                          disable=disable_progress,
                          desc="nproc="+str(nproc)):
      results.append(result)
del pool

# Sort CODFile objects depending on how they failed
failed_codfiles = [codfile for codfile in results if codfile.failed]
failed_messages = [codfile.fail_message for codfile in results if codfile.failed]
kept_codfiles = [codfile for codfile in results if not codfile.failed]
print("Failed:",len(failed_codfiles))
print("Kept:",len(kept_codfiles))


fail_counter = {message:0 for message in results[0].failure_messages.values()}
for message in failed_messages:
  fail_counter[message]+=1
  
print("\nReasons for failure:")
fail_check = 0
for key,value in fail_counter.items():
  if value>0:
    print("  ",key.ljust(60," "),":",value)
    fail_check+=value
#print("Fail check:",fail_check)


print("Files passing filters:")
for cod_file in kept_codfiles:
  print(cod_file.file)



