package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/biogo/biogo/align"
	"github.com/biogo/biogo/align/matrix"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
)

type Data struct {
	ID string
	//original
	Seq    string
	Dssp   string
	Stride string
	Kaksi  string
	Pross  string
	//processed
	Dssp3   string
	Stride3 string
	Kaksi3  string
	Pross3  string
	// consensus 2
	DsspStride3  string
	DsspKaksi3   string
	DsspPross3   string
	StrideKaksi3 string
	StridePross3 string
	KaksiPross3  string
	// consensus 3
	DsspStrideKaksi3  string
	DsspStridePross3  string
	DsspKaksiPross3   string
	StrideKaksiPross3 string
	// consensus 4
	All3 string
}

func consensus2(ss0, ss1 string) string {
	cons := ""
	for i := 0; i < len(ss0); i++ {
		if ss0[i] == ss1[i] {
			cons += string(ss0[i])
		} else {
			cons += "?"
		}
	}
	return cons
}

func consensus3(ss0, ss1, ss2 string) string {
	cons := ""
	for i := 0; i < len(ss0); i++ {
		if ss0[i] == ss1[i] && ss0[i] == ss2[i] {
			cons += string(ss0[i])
		} else {
			cons += "?"
		}
	}
	return cons
}

func consensus4(ss0, ss1, ss2, ss3 string) string {
	cons := ""
	for i := 0; i < len(ss0); i++ {
		if ss0[i] == ss1[i] && ss0[i] == ss2[i] && ss0[i] == ss3[i] {
			cons += string(ss0[i])
		} else {
			cons += "?"
		}
	}
	return cons
}

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

func dssp3(dsspss string) string {
	dssp3 := ""
	for _, v := range dsspss {
		switch v {
		case '\n':
			continue
		case 'E':
			dssp3 += "|"
		case 'G':
			dssp3 += "*"
		case 'H':
			dssp3 += "*"
		case 'I':
			dssp3 += "*"
		case 'C':
			dssp3 += "-"
		case 'B':
			dssp3 += "-"
		case 'T':
			dssp3 += "-"
		case 'S':
			dssp3 += "-"
		default:
			dssp3 += "?"
		}
	}
	return dssp3
}

func stride3(stridess string) string { return dssp3(stridess) }

func kaksi3(kaksiss string) string {
	kaksi3 := ""
	for _, v := range kaksiss {
		switch v {
		case '\n':
			continue
		case 'b':
			kaksi3 += "|"
		case 'H':
			kaksi3 += "*"
		case '.':
			kaksi3 += "-"
		default:
			kaksi3 += "?"
		}
	}
	return kaksi3
}

func pross3(prossss string) string {
	pross3 := ""
	for _, v := range prossss {
		switch v {
		case '\n':
			continue
		case 'E':
			pross3 += "|"
		case 'H':
			pross3 += "*"
		case 'T':
			pross3 += "-"
		case 'C':
			pross3 += "-"
		case 'P':
			pross3 += "-"
		default:
			pross3 += "?"
		}
	}
	return pross3
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
	// To debug
	// fmt.Println(seq_aa)
	// fmt.Println(dssp_aa)
	// fmt.Println(stride_aa)
	// fmt.Println(pross_aa)
	// fmt.Println(kaksi_aa)

	// Monta os alinhamentos
	// aa, dsspaa, dsspss := alinha(seq_aa, dssp_aa, dssp_ss)
	aa, _, dsspss := alinha(seq_aa, dssp_aa, dssp_ss)
	// _, strideaa, stridess := alinha(seq_aa, stride_aa, stride_ss)
	_, _, stridess := alinha(seq_aa, stride_aa, stride_ss)
	// _, prossaa, prossss := alinha(seq_aa, pross_aa, pross_ss)
	_, _, prossss := alinha(seq_aa, pross_aa, pross_ss)
	// _, kaksiaa, kaksiss := alinha(seq_aa, kaksi_aa, kaksi_ss)
	_, _, kaksiss := alinha(seq_aa, kaksi_aa, kaksi_ss)

	// fmt.Printf("PDBSEQUENCE %s\n", aa)
	// // fmt.Printf("DSSPSEQ     %s\n", dsspaa)
	// fmt.Printf("DSSPSECSTR  %s\n", dsspss)
	// // fmt.Printf("STRDSEQ     %s\n", strideaa)
	// fmt.Printf("STRDSECSTR  %s\n", stridess)
	// // fmt.Printf("PROSSEQ     %s\n", prossaa)
	// fmt.Printf("PROSSECSTR  %s\n", prossss)
	// // fmt.Printf("KAKSSEQ     %s\n", kaksiaa)
	// fmt.Printf("KAKSSECSTR  %s\n", kaksiss)

	// Escreve um arquivo com os dados de ss e gaps em um arquivo JSON
	id := filepath.Base(*fastaFName)
	id = strings.TrimSuffix(id, ".fa")
	id = strings.TrimSuffix(id, ".fasta")

	data := Data{
		ID:                id,
		Seq:               aa,
		Dssp:              dsspss,
		Stride:            stridess,
		Kaksi:             kaksiss,
		Pross:             prossss,
		Dssp3:             dssp3(dsspss),
		Stride3:           stride3(stridess),
		Kaksi3:            kaksi3(kaksiss),
		Pross3:            pross3(prossss),
		DsspStride3:       consensus2(dssp3(dsspss), stride3(stridess)),
		DsspKaksi3:        consensus2(dssp3(dsspss), kaksi3(kaksiss)),
		DsspPross3:        consensus2(dssp3(dsspss), pross3(prossss)),
		StrideKaksi3:      consensus2(stride3(stridess), kaksi3(kaksiss)),
		StridePross3:      consensus2(stride3(stridess), pross3(prossss)),
		KaksiPross3:       consensus2(kaksi3(kaksiss), pross3(prossss)),
		DsspStrideKaksi3:  consensus3(dssp3(dsspss), stride3(stridess), kaksi3(kaksiss)),
		DsspStridePross3:  consensus3(dssp3(dsspss), stride3(stridess), pross3(prossss)),
		DsspKaksiPross3:   consensus3(dssp3(dsspss), kaksi3(kaksiss), pross3(prossss)),
		StrideKaksiPross3: consensus3(stride3(stridess), kaksi3(kaksiss), pross3(prossss)),
		All3:              consensus4(dssp3(dsspss), stride3(stridess), kaksi3(kaksiss), pross3(prossss)),
	}
	b, err := json.MarshalIndent(data, "", "  ")
	if err != nil {
		fmt.Println("error:", err)
	}
	fmt.Println(string(b))
}
