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
		f_aa = align.Format(seq_aa, ss_aa, aln, '.')
		f_ss = align.Format(seq_aa, ss_ss, aln, '.')
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
			dssp3 += "_"
		case 'B':
			dssp3 += "_"
		case 'T':
			dssp3 += "_"
		case 'S':
			dssp3 += "_"
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
			kaksi3 += "_"
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
			pross3 += "_"
		case 'C':
			pross3 += "_"
		case 'P':
			pross3 += "_"
		default:
			pross3 += "?"
		}
	}
	return pross3
}

func hydro(seq, scale string) string {
	var hp string
	if scale == "rose" {
		for _, v := range seq {
			hp += rose[string(v)]
		}
	}

	if scale == "roseSpecial" {
		for _, v := range seq {
			hp += roseSpecial[string(v)]
		}
	}

	if scale == "roseSpecialCharged" {
		for _, v := range seq {
			hp += roseSpecialCharged[string(v)]
		}
	}

	return hp
}

func combine(seq, hp string) []string {
	cb := make([]string, len(seq))
	for i := 0; i < len(cb); i++ {
		cb[i] = string(seq[i]) + string(hp[i])
	}
	return cb

}

// func roseHP(seq string) string {
// 	hp := ""
// 	for _, v := range seq {
// 		switch v {
// 		case 'D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Y', 'G', 'P':
// 			hp += "p"
// 		case 'A', 'C', 'L', 'I', 'F', 'W', 'V', 'M':
// 			hp += "n"
// 		default:
// 			fmt.Println("Residue not found", v)
// 			return ""
// 		}
// 	}
// 	return hp
// }
//
// func roseHPSpecialPG(seq string) string {
// 	hp := ""
// 	for _, v := range seq {
// 		switch v {
// 		case 'D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Y':
// 			hp += "p"
// 		case 'A', 'C', 'L', 'I', 'F', 'W', 'V', 'M':
// 			hp += "n"
// 		case 'P':
// 			hp += "P"
// 		case 'G':
// 			hp += "G"
// 		default:
// 			fmt.Println("Residue not found", v)
// 			return ""
// 		}
// 	}
// 	return hp
// }

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

	hpRose := hydro(aa, "rose")
	hpRoseSpecial := hydro(aa, "roseSpecial")
	hpRoseSpecialCharged := hydro(aa, "roseSpecialCharged")

	data := Data{
		ID:                id,
		Seq:               strings.Split(aa, ""),
		Dssp:              strings.Split(dsspss, ""),
		Stride:            strings.Split(stridess, ""),
		Kaksi:             strings.Split(kaksiss, ""),
		Pross:             strings.Split(prossss, ""),
		Dssp3:             strings.Split(dssp3(dsspss), ""),
		Stride3:           strings.Split(stride3(stridess), ""),
		Kaksi3:            strings.Split(kaksi3(kaksiss), ""),
		Pross3:            strings.Split(pross3(prossss), ""),
		DsspStride3:       strings.Split(consensus2(dssp3(dsspss), stride3(stridess)), ""),
		DsspKaksi3:        strings.Split(consensus2(dssp3(dsspss), kaksi3(kaksiss)), ""),
		DsspPross3:        strings.Split(consensus2(dssp3(dsspss), pross3(prossss)), ""),
		StrideKaksi3:      strings.Split(consensus2(stride3(stridess), kaksi3(kaksiss)), ""),
		StridePross3:      strings.Split(consensus2(stride3(stridess), pross3(prossss)), ""),
		KaksiPross3:       strings.Split(consensus2(kaksi3(kaksiss), pross3(prossss)), ""),
		DsspStrideKaksi3:  strings.Split(consensus3(dssp3(dsspss), stride3(stridess), kaksi3(kaksiss)), ""),
		DsspStridePross3:  strings.Split(consensus3(dssp3(dsspss), stride3(stridess), pross3(prossss)), ""),
		DsspKaksiPross3:   strings.Split(consensus3(dssp3(dsspss), kaksi3(kaksiss), pross3(prossss)), ""),
		StrideKaksiPross3: strings.Split(consensus3(stride3(stridess), kaksi3(kaksiss), pross3(prossss)), ""),
		All3:              strings.Split(consensus4(dssp3(dsspss), stride3(stridess), kaksi3(kaksiss), pross3(prossss)), ""),
		//rose hydrophobicity
		SeqHPRose:               combine(aa, hpRose),
		DsspHPRose:              combine(dsspss, hpRose),
		StrideHPRose:            combine(stridess, hpRose),
		KaksiHPRose:             combine(kaksiss, hpRose),
		ProssHPRose:             combine(prossss, hpRose),
		Dssp3HPRose:             combine(dssp3(dsspss), hpRose),
		Stride3HPRose:           combine(stride3(stridess), hpRose),
		Kaksi3HPRose:            combine(kaksi3(kaksiss), hpRose),
		Pross3HPRose:            combine(pross3(prossss), hpRose),
		DsspStride3HPRose:       combine(consensus2(dssp3(dsspss), stride3(stridess)), hpRose),
		DsspKaksi3HPRose:        combine(consensus2(dssp3(dsspss), kaksi3(kaksiss)), hpRose),
		DsspPross3HPRose:        combine(consensus2(dssp3(dsspss), pross3(prossss)), hpRose),
		StrideKaksi3HPRose:      combine(consensus2(stride3(stridess), kaksi3(kaksiss)), hpRose),
		StridePross3HPRose:      combine(consensus2(stride3(stridess), pross3(prossss)), hpRose),
		KaksiPross3HPRose:       combine(consensus2(kaksi3(kaksiss), pross3(prossss)), hpRose),
		DsspStrideKaksi3HPRose:  combine(consensus3(dssp3(dsspss), stride3(stridess), kaksi3(kaksiss)), hpRose),
		DsspStridePross3HPRose:  combine(consensus3(dssp3(dsspss), stride3(stridess), pross3(prossss)), hpRose),
		DsspKaksiPross3HPRose:   combine(consensus3(dssp3(dsspss), kaksi3(kaksiss), pross3(prossss)), hpRose),
		StrideKaksiPross3HPRose: combine(consensus3(stride3(stridess), kaksi3(kaksiss), pross3(prossss)), hpRose),
		All3HPRose:              combine(consensus4(dssp3(dsspss), stride3(stridess), kaksi3(kaksiss), pross3(prossss)), hpRose),

		//rose hydrophobicity
		SeqHPRoseSpecial:               combine(aa, hpRoseSpecial),
		DsspHPRoseSpecial:              combine(dsspss, hpRoseSpecial),
		StrideHPRoseSpecial:            combine(stridess, hpRoseSpecial),
		KaksiHPRoseSpecial:             combine(kaksiss, hpRoseSpecial),
		ProssHPRoseSpecial:             combine(prossss, hpRoseSpecial),
		Dssp3HPRoseSpecial:             combine(dssp3(dsspss), hpRoseSpecial),
		Stride3HPRoseSpecial:           combine(stride3(stridess), hpRoseSpecial),
		Kaksi3HPRoseSpecial:            combine(kaksi3(kaksiss), hpRoseSpecial),
		Pross3HPRoseSpecial:            combine(pross3(prossss), hpRoseSpecial),
		DsspStride3HPRoseSpecial:       combine(consensus2(dssp3(dsspss), stride3(stridess)), hpRoseSpecial),
		DsspKaksi3HPRoseSpecial:        combine(consensus2(dssp3(dsspss), kaksi3(kaksiss)), hpRoseSpecial),
		DsspPross3HPRoseSpecial:        combine(consensus2(dssp3(dsspss), pross3(prossss)), hpRoseSpecial),
		StrideKaksi3HPRoseSpecial:      combine(consensus2(stride3(stridess), kaksi3(kaksiss)), hpRoseSpecial),
		StridePross3HPRoseSpecial:      combine(consensus2(stride3(stridess), pross3(prossss)), hpRoseSpecial),
		KaksiPross3HPRoseSpecial:       combine(consensus2(kaksi3(kaksiss), pross3(prossss)), hpRoseSpecial),
		DsspStrideKaksi3HPRoseSpecial:  combine(consensus3(dssp3(dsspss), stride3(stridess), kaksi3(kaksiss)), hpRoseSpecial),
		DsspStridePross3HPRoseSpecial:  combine(consensus3(dssp3(dsspss), stride3(stridess), pross3(prossss)), hpRoseSpecial),
		DsspKaksiPross3HPRoseSpecial:   combine(consensus3(dssp3(dsspss), kaksi3(kaksiss), pross3(prossss)), hpRoseSpecial),
		StrideKaksiPross3HPRoseSpecial: combine(consensus3(stride3(stridess), kaksi3(kaksiss), pross3(prossss)), hpRoseSpecial),
		All3HPRoseSpecial:              combine(consensus4(dssp3(dsspss), stride3(stridess), kaksi3(kaksiss), pross3(prossss)), hpRoseSpecial),

		//rose hydrophobicity
		SeqHPRoseSpecialCharged:               combine(aa, hpRoseSpecialCharged),
		DsspHPRoseSpecialCharged:              combine(dsspss, hpRoseSpecialCharged),
		StrideHPRoseSpecialCharged:            combine(stridess, hpRoseSpecialCharged),
		KaksiHPRoseSpecialCharged:             combine(kaksiss, hpRoseSpecialCharged),
		ProssHPRoseSpecialCharged:             combine(prossss, hpRoseSpecialCharged),
		Dssp3HPRoseSpecialCharged:             combine(dssp3(dsspss), hpRoseSpecialCharged),
		Stride3HPRoseSpecialCharged:           combine(stride3(stridess), hpRoseSpecialCharged),
		Kaksi3HPRoseSpecialCharged:            combine(kaksi3(kaksiss), hpRoseSpecialCharged),
		Pross3HPRoseSpecialCharged:            combine(pross3(prossss), hpRoseSpecialCharged),
		DsspStride3HPRoseSpecialCharged:       combine(consensus2(dssp3(dsspss), stride3(stridess)), hpRoseSpecialCharged),
		DsspKaksi3HPRoseSpecialCharged:        combine(consensus2(dssp3(dsspss), kaksi3(kaksiss)), hpRoseSpecialCharged),
		DsspPross3HPRoseSpecialCharged:        combine(consensus2(dssp3(dsspss), pross3(prossss)), hpRoseSpecialCharged),
		StrideKaksi3HPRoseSpecialCharged:      combine(consensus2(stride3(stridess), kaksi3(kaksiss)), hpRoseSpecialCharged),
		StridePross3HPRoseSpecialCharged:      combine(consensus2(stride3(stridess), pross3(prossss)), hpRoseSpecialCharged),
		KaksiPross3HPRoseSpecialCharged:       combine(consensus2(kaksi3(kaksiss), pross3(prossss)), hpRoseSpecialCharged),
		DsspStrideKaksi3HPRoseSpecialCharged:  combine(consensus3(dssp3(dsspss), stride3(stridess), kaksi3(kaksiss)), hpRoseSpecialCharged),
		DsspStridePross3HPRoseSpecialCharged:  combine(consensus3(dssp3(dsspss), stride3(stridess), pross3(prossss)), hpRoseSpecialCharged),
		DsspKaksiPross3HPRoseSpecialCharged:   combine(consensus3(dssp3(dsspss), kaksi3(kaksiss), pross3(prossss)), hpRoseSpecialCharged),
		StrideKaksiPross3HPRoseSpecialCharged: combine(consensus3(stride3(stridess), kaksi3(kaksiss), pross3(prossss)), hpRoseSpecialCharged),
		All3HPRoseSpecialCharged:              combine(consensus4(dssp3(dsspss), stride3(stridess), kaksi3(kaksiss), pross3(prossss)), hpRoseSpecialCharged),
	}
	b, err := json.MarshalIndent(data, "", "  ")
	if err != nil {
		fmt.Println("error:", err)
	}
	fmt.Println(string(b))
}
