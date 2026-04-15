from pyfaidx import Fasta
import editdistance
import re


def define_umi_simple(file, output, filename):
	"""
	Use this function if the UMI is the first part of the read. e.g.
	<UMI> <UMI boundary> <Read>
	"""
	fasta_records = Fasta(file)
	filename = filename = '.'.join(filename.split('.')[:-1])
	
	with open(output+filename+'_umi.fasta', 'w') as f:
	
		for record_id, sequence in fasta_records.items():
			umi = fasta_records[record_id][:12]

			f.write('>' + record_id + '_' + str(umi) + '\n')
			f.write(str(sequence)[12:] + '\n')

			
def define_umi(file, output, filename):
	"""
	Use this function if the adapter is kept during C3POa post-processing. Change adapter string to match. e.g.:
		<adapter> <UMI> <UMI boundary> <Read>
	"""
	adapter = 'TGTATAAGAGACAG'

	fasta_records = Fasta(file)
	filename = filename.split('.')[0]
	
	with open(output + filename + '_umi.fasta', 'w') as f:
		
		for record_id, sequence in fasta_records.items():
			read_start = fasta_records[record_id][:73]  # Depends on C3POa Adapter (was 33 with no de-multiplexing/73 with demux)
			start_index = 32  # How many bases in the start of the adapter sequence // 32 with demux, 7 with no demux
			slice_length = 14  # Length of the adapter sequence // 14 for both
			aligned_segment = read_start[start_index:start_index + slice_length]  # Depends on C3POa Adapter
			if editdistance.eval(adapter, str(aligned_segment)) <= 4:  # 2 with no demux // 4 with demux
				tmp_umi = re.search('ACAG' + '.{12}', str(read_start))
				if tmp_umi is None:
					continue
				umi = tmp_umi.group(0)[4:]
				
				f.write('>' + record_id + '_' + umi + '\n')
				f.write(str(sequence) + '\n')  # Depends on C3POa Adapter (was 33 with no de-multiplexing/73 with demux)
				