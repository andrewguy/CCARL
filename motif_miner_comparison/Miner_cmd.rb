require 'json'
require './khashimo/Require'
require_relative 'lib/CMTreeMiner'

minsup = ARGV[0]
alpha = ARGV[1]
glycan_file = ARGV[2]

glycan_data = File.read(glycan_file)

cmt = CMTreeMiner.new(minsup.to_i, K_Entry_arr.new(glycan_data), alpha.to_f)
cmt.mining("closed")

for subtree in cmt.getSubtrees do
    puts subtree.getKCF
end
