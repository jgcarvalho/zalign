package main

import (
	"fmt"
	"os"

	"github.com/biogo/biogo/align"
	"github.com/biogo/biogo/align/matrix"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
	//	"github.com/biogo/biogo/alphabet"
	"flag"

	//"github.com/biogo/biogo/io/seqio/fasta"
	//	"github.com/biogo/biogo/seq/linear"
	//	"reflect"

	//	"pepplanes"
	// "github.com/jgcarvalho/zutils/ssparser"
)

func getFasta(fn string) (seq.Sequence, error) {
	fasta_file, err := os.Open(fn)
	if err != nil {
		fmt.Println("Erro ao ler o arquivo", err)
	}
	var s []alphabet.Letter
	t := linear.NewSeq("", s, alphabet.Protein)
	reader := fasta.NewReader(fasta_file, t)
	seq, _ := reader.Read()
	//fmt.Println("Read -> ", seq.Alphabet())
	return seq, nil
}

func getAttr(fn string) (seq.Sequence, seq.Sequence, error) {
	fasta_file, err := os.Open(fn)
	if err != nil {
		fmt.Println("Erro ao ler o arquivo", err)
	}
	var s []alphabet.Letter
	t := linear.NewSeq("", s, alphabet.Protein)
	reader := fasta.NewReader(fasta_file, t)
	seq, _ := reader.Read()
	// if reader.Next() != nil {
	// 	fmt.Println("PAU NO NEXT")
	// }
	attr, _ := reader.Read()
	//fmt.Println("Read -> ", seq.Alphabet())
	return seq, attr, nil
}

func alinha(seq_aa, ss_aa, ss_ss seq.Sequence) (string, string, string) {
	nw := align.NWAffine{
		Matrix:  matrix.MATCH,
		GapOpen: -1,
	}

	aln, err := nw.Align(seq_aa, ss_aa)
	var f_aa, f_ss [2]alphabet.Slice
	if err == nil {
		f_aa = align.Format(seq_aa, ss_aa, aln, '_')
		f_ss = align.Format(seq_aa, ss_ss, aln, '_')
	} else {
		fmt.Println("O ERRO E:", err)
	}
	return fmt.Sprint(f_aa[0]), fmt.Sprint(f_aa[1]), fmt.Sprint(f_ss[1])
}

func main() {

	fastaFName := flag.String("fa", "", "fasta file")
	// // pdbFName := flag.String("pdb", "", "pdb file")
	dsspFName := flag.String("dssp", "", "dssp file")
	strideFName := flag.String("stride", "", "stride file")
	prossFName := flag.String("pross", "", "pross file")
	kaksiFName := flag.String("kaksi", "", "kaksi file")
	// // doKaksi := flag.Bool("kaksi", false, "kaksi")
	// // doBba := flag.Bool("bba", false, "bba")
	// // classBba := flag.Int("class_bba", 0, "class bba")
	// // classAa := flag.Int("class_aa", 0, "class aa")

	flag.Parse()

	var seq_aa,
		dssp_aa,
		dssp_ss,
		stride_aa,
		stride_ss,
		pross_aa,
		pross_ss,
		kaksi_aa,
		kaksi_ss seq.Sequence
	var err error

	// Le a sequencia fasta
	if *fastaFName == "" {
		fmt.Println("É necessario o arquivo fasta")
		return
	} else {
		// seq_aa, err = getFasta("/home/jgcarvalho/sync/ZZpred/pisces/fasta/1A1XA.fa")
		seq_aa, err = getFasta(*fastaFName)
		if err != nil {
			fmt.Println("Erro no processamento do Fasta", err)
		}
	}

	if *dsspFName == "" {
		fmt.Println("É necessario o arquivo DSSP")
		return
	} else {
		// dssp_aa, dssp_ss, err = getAttr("/home/jgcarvalho/sync/ZZpred/pisces/dssp_chain_fasta/1A1XA_dssp.fasta")
		dssp_aa, dssp_ss, err = getAttr(*dsspFName)
		if err != nil {
			fmt.Println("Erro no processamento do DSSP", err)
		}
	}

	if *strideFName == "" {
		fmt.Println("É necessario o arquivo STRIDE")
		return
	} else {
		// stride_aa, stride_ss, err = getAttr("/home/jgcarvalho/sync/ZZpred/pisces/stride_chain_fasta/1A1XA_stride.fasta")
		stride_aa, stride_ss, err = getAttr(*strideFName)
		if err != nil {
			fmt.Println("Erro no processamento do STRIDE", err)
		}
	}

	if *prossFName == "" {
		fmt.Println("É necessario o arquivo PROSS")
		return
	} else {
		// stride_aa, stride_ss, err = getAttr("/home/jgcarvalho/sync/ZZpred/pisces/stride_chain_fasta/1A1XA_stride.fasta")
		pross_aa, pross_ss, err = getAttr(*prossFName)
		if err != nil {
			fmt.Println("Erro no processamento do STRIDE", err)
		}
	}

	if *kaksiFName == "" {
		fmt.Println("É necessario o arquivo KAKSI")
		return
	} else {
		// stride_aa, stride_ss, err = getAttr("/home/jgcarvalho/sync/ZZpred/pisces/stride_chain_fasta/1A1XA_stride.fasta")
		kaksi_aa, kaksi_ss, err = getAttr(*kaksiFName)
		if err != nil {
			fmt.Println("Erro no processamento do KAKSI", err)
		}
	}

	// Monta os alinhamentos
	// aa, dsspaa, dsspss := alinha(seq_aa, dssp_aa, dssp_ss)
	aa, _, dsspss := alinha(seq_aa, dssp_aa, dssp_ss)
	// _, strideaa, stridess := alinha(seq_aa, stride_aa, stride_ss)
	_, _, stridess := alinha(seq_aa, stride_aa, stride_ss)
	// _, prossaa, prossss := alinha(seq_aa, pross_aa, pross_ss)
	_, _, prossss := alinha(seq_aa, pross_aa, pross_ss)
	// _, kaksiaa, kaksiss := alinha(seq_aa, kaksi_aa, kaksi_ss)
	_, _, kaksiss := alinha(seq_aa, kaksi_aa, kaksi_ss)
	fmt.Printf("PDBSEQUENCE %s\n", aa)
	// fmt.Printf("DSSPSEQ     %s\n", dsspaa)
	fmt.Printf("DSSPSECSTR  %s\n", dsspss)
	// fmt.Printf("STRDSEQ     %s\n", strideaa)
	fmt.Printf("STRDSECSTR  %s\n", stridess)
	// fmt.Printf("PROSSEQ     %s\n", prossaa)
	fmt.Printf("PROSSECSTR  %s\n", prossss)
	// fmt.Printf("KAKSSEQ     %s\n", kaksiaa)
	fmt.Printf("KAKSSECSTR  %s\n", kaksiss)

	// Escreve um arquivo com os dados de ss e gaps em um arquivo JSON
}
