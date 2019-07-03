from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

import googletrans
from translate import Translator
#from googletrans import Translator
#translator = Translator()


def normalized_distance(word1, word2):
	if word1 == word2:
		dis = 1
	else:
		#dis = pairwise2.align.globalxx(word1, word2, score_only=True, one_alignment_only=True)
		dis = pairwise2.align.globalmx(list(word1), list(word2), 1, -0.5, gap_char=['-'], score_only=True, one_alignment_only=True)
		dis = dis / max( len(word1), len(word2) )
	return dis

def max_pairwise_align(w1, w2):

	#print(translator.translate('house', dest='it').text.lower())
	#print(googletrans.LANGUAGES)

	# Manual deletions to avoid core dump
	languages = googletrans.LANGUAGES
	del languages['ar']

	# Build dict pairing LANGUAGE/matrix index
	index = {}
	index_lang = {}
	i = 0
	for (item, lang) in languages.items():
		index[item] = i
		index_lang[i] = item
		i = i + 1
	
	# Build MxM result matrix
	W = len( languages.items() )
	Dis = [[0 for x in range(W)] for y in range(W)] 

	# Run translation
	words = ["was", "many", "make", "long"]
	word = "car"
	src = "en"
	Tran = index.copy()
	Tran[src] = word
	for (item, lang) in languages.items():
		if item != src:
			#Tran[item] = translator.translate(word, src=src, dest=item).text.lower()
			translator = Translator(from_lang=src, to_lang=item)
			Tran[item] = translator.translate(word) 

	#print(Tran)
	# Calculate distances
	for i in range(W-1):
		for j in range(i+1, W):
			lang1 = index_lang[i]
			lang2 = index_lang[j]
			word1 = Tran[lang1]
			word2 = Tran[lang2]
			print(i, j, lang1, lang2, word1, word2)
			try:
				distance = normalized_distance(word1, word2)
			except Excetion:
				distance = 0
			Dis[j][i] = Dis[j][i] + distance
			
	Dis2 = []
	for x in range(W):
		Dis2.append( Dis[x][0:x] )

	dm = DistanceMatrix(names=languages.keys(), matrix=Dim2)
	draw_tree(dm)


def draw_tree(distance_matrix, model='upgma'):
	constructor = DistanceTreeConstructor()
	if model == 'upgma':
		tree = constructor.upgma(dm)
	elif model == 'nj':
		tree = constructor.nj(dm)
	else:
		print("Error model")
		exit(1)
	print(tree)
	Phylo.draw(tree)

max_pairwise_align("a", "aa")
#draw_tree(dm)