#!/usr/bin/env ruby

require 'bio'

########################################
#This program parses out the core regions from a mauve alignment
#It won't work on machines that didn't run the Mauve job without altering
#the file locations in the xmfa file.
#
#Usage: ruby extract_core_regions.rb xmfa_file backbone_file
#
########################################

xmfa_file = ARGV[0]
backbone_file = ARGV[1]
core_regions_file = "core_regions.tsv"

if File.exists?(xmfa_file) && File.file?(xmfa_file) 
then
  xmfa_FH = File.open(xmfa_file, "r")
else
  puts "Problem opening #{xmfa_file}"
  exit 1
end

#Get which sequences are which (i.e. seq1 refers to which fasta file?)
sequence_dictionary = Array.new

xmfa_FH.each do |line|
  if /\#Backbone/.match(line)
    break
  elsif my_match = /\#Sequence(\d+)File\s+(.+)/.match(line)
    update_num = my_match[1].to_i-1
    sequence_dictionary << [update_num, my_match[2]]
  end
end

#Now, create our modified backbone file, with only core sequences
                
if File.exists?(backbone_file) && File.file?(backbone_file)
  core_regions = %x[grep -v '[[:space:]]0[[:space:]]' #{backbone_file}]
  out_FH = File.open(core_regions_file, "w")
  out_FH.write(core_regions)
  out_FH.close
else 
  puts "Problem opening #{backbone_file}"
  exit 1
end


#Iterate through the newly created core regions file, and add all the positions to the matrix
core_position_matrix = Array.new

File.open(core_regions_file, "r").each do |record|
  core_position_matrix << record.split("\t")
end

#puts core_position_matrix.size

#Finally we can open the sequence files and get the positions, this is the meat of the program
sequence_dictionary.each do |seq_file|
  
#Need the index of our current sequence in the matrix
  left = 0
  right = 0
  
  #Get the _leftend & _rightend column positions for each genome
  core_position_matrix[0].each do |seq_end|

    if seq_end == "seq"+seq_file[0].to_s+"_leftend"
      left = core_position_matrix[0].index(seq_end).to_i
    end    

  end
  
  right = left + 1;

  #The backbone file gives positions for a "concatenated" file, so we need to concatenate all the contigs of each genome
  seq_file[1].chomp!
  if File.exists?(seq_file[1]) && File.file?(seq_file[1])
    #There is a problem sometimes where the "auto" function doesn't recognize RAST genbanks.  You can force it by uncommenting the next line, and commenting out the line after:
    #seq_FH = Bio::FlatFile.open(Bio::GenBank, seq_file[1])
    seq_FH = Bio::FlatFile.auto(seq_file[1])
    full_seq = ""
    seq_FH.each do |seq_entry|
      full_seq << seq_entry.seq
    end
  else
    puts "Problem opening "+seq_file[1]
    exit 1
  end
  
  bio_full_seq = Bio::Sequence::NA.new(full_seq)
  core_seq = ""
  
  core_position_matrix.drop(1).each do |core_position|

     #Logic to handle rev-comp sequences, and normal sequences
    if core_position[left].to_i < 0 || core_position[right].to_i < 0
      core_seq << Bio::Sequence::NA.new(bio_full_seq.subseq(core_position[left].to_i.abs, core_position[right].to_i.abs)).complement
    else
      #Some error handling junk to make sure we are actually getting sequence matching the positions
      #TODO: this is a terrible spot for this check
      if bio_full_seq.subseq(core_position[left].to_i.abs, core_position[right].to_i.abs).nil?
        puts core_position[left]+" "+core_position[right]+" resolves to nil"
      else
        core_seq << bio_full_seq.subseq(core_position[left].to_i.abs, core_position[right].to_i.abs)
      end
    end
  end
  puts ">"+seq_file[1]
  puts core_seq
#  puts "Sequence is identified by mauve as seq"+seq_file[0].to_s+" the left side is: "+core_position_matrix[0][left]+" and the right is "+core_position_matrix[0][right]
#  puts "Sequence "+seq_file[1].to_s+" has a core region size of "+core_seq.size.to_s+" And a total size of "+full_seq.size.to_s

end


