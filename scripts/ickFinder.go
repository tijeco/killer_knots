package main

import (
	"bufio"
	"bytes"
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"regexp"
	"strconv"
	"strings"
	"time"
)

var p = fmt.Println

var ickDB, queryPep, outDir, signalPath string
var runSilix, help bool
var iterate int

func init() {

	flag.StringVar(&ickDB, "db", "", "fasta file with reference ICK")
	flag.StringVar(&queryPep, "pep", "", "fasta file with query peptides")
	flag.StringVar(&outDir, "out", "", "output directory")
	flag.StringVar(&outDir, "o", "", "output directory (shorthand)")
	flag.StringVar(&signalPath, "signalp", "", "full path to signalP binary")
	flag.BoolVar(&runSilix, "silix", false, "run silix on final output to cluster results")
	flag.IntVar(&iterate, "n", 1, "Number of times to iterate")
	flag.BoolVar(&help, "help", false, "print usage")
	flag.BoolVar(&help, "h", false, "print usage (shorthand)")
	flag.Parse()

	someEmpty := ickDB == "" || queryPep == "" || outDir == ""
	// p(someEmpty)
	if someEmpty {
		flag.PrintDefaults()
		os.Exit(1)
	} else {
		outerr := os.MkdirAll(outDir, os.ModePerm)
		if outerr != nil {
			log.Fatal(outerr)
		}
		for _, file := range []string{ickDB, queryPep} {
			if fileExists(file) == false {
				log.Fatalln(file, "does not exist")
			}
		}

	}
	// p("checking if signalP has been requested")
	// log.Fatalln(signalPath == "")

}

type fasta struct {
	id   string
	desc string
	seq  string
}

func timeStatus(status string) {
	p(time.Now().Format("15:04:05 Mon Jan-02-2006"), "    ", status)
	p()
}

func logFatalErr(message string, err error) {
	if err != nil {
		log.Fatal(err)
	}
}

func checkDependencies(dependencies []string) (good2go bool) {
	for _, dependency := range dependencies {
		good2go = willRun(dependency)
	}
	return
}

func splitLines(s string) []string {
	var lines []string
	sc := bufio.NewScanner(strings.NewReader(s))
	for sc.Scan() {
		lines = append(lines, sc.Text())
	}
	return lines
}

func buildFasta(header string, seq bytes.Buffer) (record fasta) {
	fields := strings.SplitN(header, " ", 2)

	if len(fields) > 1 {
		record.id = fields[0]
		record.desc = fields[1]
	} else {
		record.id = fields[0]
		record.desc = ""
	}

	record.seq = strings.ToUpper(seq.String())

	return record
}

func parseFasta(fastaFh io.Reader) chan fasta {

	outputChannel := make(chan fasta)

	scanner := bufio.NewScanner(fastaFh)
	// scanner.Split(bufio.ScanLines)
	header := ""
	var seq bytes.Buffer

	go func() {
		// Loop over the letters in inputString
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if len(line) == 0 {
				continue
			}

			// line := scanner.Text()

			if line[0] == '>' {
				// If we stored a previous identifier, get the DNA string and map to the
				// identifier and clear the string
				if header != "" {
					// outputChannel <- buildFasta(header, seq.String())
					outputChannel <- buildFasta(header, seq)
					// fmt.Println(record.id, len(record.seq))
					header = ""
					seq.Reset()
				}

				// Standard FASTA identifiers look like: ">id desc"
				header = line[1:]
			} else {
				// Append here since multi-line DNA strings are possible
				seq.WriteString(line)
			}

		}

		outputChannel <- buildFasta(header, seq)

		// Close the output channel, so anything that loops over it
		// will know that it is finished.
		close(outputChannel)
	}()

	return outputChannel
}

// try using it to prevent further errors.
func fileExists(filename string) bool {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}

func willRun(command string) (ran bool) {
	commandPath, err := exec.LookPath(command)
	timeStatus("Checking if " + command + " is installed")
	if err != nil {
		log.Fatal(command + " not found in $PATH, exiting now")
	} else {
		timeStatus(command + " is available at: " + commandPath)
		ran = true
	}
	return
}

func blastp(query, target, outDir, outName string, onlyOne bool, evalue float64) (blastTargets []string) {
	var blastP string
	blastOutFile := outDir + "/" + outName + ".txt"
	if fileExists(blastOutFile) {
		timeStatus(blastOutFile + " exists, opening it now")
		blastPstdout, blastPerr := ioutil.ReadFile(blastOutFile)
		if blastPerr != nil {
			log.Fatal(blastPerr)
		} else {
			blastTargets = targetFromBlast(string(blastPstdout))
		}
	} else {
		// if willRun("makeblastdb") && willRun("blastp") {
		// makeblastdb -in spiderICK.fasta  -dbtype prot
		makeBlastDB := "makeblastdb -in " + target + " -dbtype prot"
		makeBlastDBcmd := exec.Command("sh", "-c", makeBlastDB)
		_, makeBlastDBerr := makeBlastDBcmd.Output()
		timeStatus(makeBlastDB)
		// p("running:", makeBlastDB)
		if makeBlastDBerr != nil {
			p(":(")
			log.Fatal(makeBlastDBerr)
		} else {
			// p("WHOOP!", makeBlastDBstdout)
			// -query Stegotoxin.protein.fa -db spiderICK.fasta -outfmt 6
			if onlyOne {
				blastP = "blastp -query " + query + " -db " + target + " -outfmt 6 -max_target_seqs 1 -evalue " + fmt.Sprintf("%f", evalue)
			} else {
				blastP = "blastp -query " + query + " -db " + target + " -outfmt 6 -evalue " + fmt.Sprintf("%f", evalue)
			}

			blastPcmd := exec.Command("sh", "-c", blastP)
			blastPstdout, blastPerr := blastPcmd.Output()
			if blastPerr != nil {
				log.Fatal(blastPerr)
			} else {
				timeStatus(blastP)
				// p("running", blastP)
				// p("WHOOP!", blastPstdout)
				blastTargets = targetFromBlast(string(blastPstdout))

				blastOutErr := ioutil.WriteFile(blastOutFile, blastPstdout, 0644)
				if blastOutErr != nil {
					log.Fatal(blastOutErr)
				}
			}

		}
	}

	// }
	return
}

func targetFromBlast(blastout string) (targets []string) {
	// p(strings.Split(blastout, "\n"))
	for _, row := range strings.Split(blastout, "\n") {
		columns := strings.Fields(row)
		// p(row)
		// p(len(columns))
		if len(columns) == 12 {
			targets = append(targets, columns[0])
		}
		// targets = append(targets, strings.Fields(row)[0])
	}
	return
}

func mafft(inFasta string) (alignment string) {
	// if willRun("mafft") {
	mafftStr := "mafft --localpair --maxiterate 1000 " + inFasta
	mafftCmd := exec.Command("sh", "-c", mafftStr)
	timeStatus(mafftStr)
	// p("running", mafftStr)
	mafftStdout, mafftErr := mafftCmd.Output()
	if mafftErr != nil {
		log.Fatal(mafftErr)
	} else {

		// p("WHOOP!", blastPstdout)
		alignment = string(mafftStdout)
	}
	// }
	return
}

func hmmer(query, target, outDir string, n int) (hmmerTargets []string) {
	// if willRun("hmmbuild") && willRun("hmmsearch") {
	var msa []byte
	var msaOutErr error
	msaOutFile := outDir + "/ickDB_" + strconv.Itoa(n) + ".aln"
	if fileExists(msaOutFile) == false {
		msa = []byte(mafft(target))
		msaOutErr = ioutil.WriteFile(msaOutFile, msa, 0644)
	} else {
		timeStatus(msaOutFile + " already exists, will use it")
	}

	if msaOutErr != nil {
		log.Fatal(msaOutErr)
	} else {
		hmmOutFile := msaOutFile + ".hmm"
		hmmbuild := "hmmbuild " + hmmOutFile + " " + msaOutFile
		hmmbuildCmd := exec.Command("sh", "-c", hmmbuild)
		timeStatus(hmmbuild)
		// p("running", hmmbuild)
		_, hmmbuildErr := hmmbuildCmd.Output()
		if hmmbuildErr != nil {
			log.Fatal(hmmbuildErr)
		} else {
			hmmsearch := "hmmsearch --notextw " + hmmOutFile + " " + query
			hmmsearchCmd := exec.Command("sh", "-c", hmmsearch)
			timeStatus(hmmsearch)
			// p("running", hmmsearch)
			hmmsearchOut, hmmsearchErr := hmmsearchCmd.Output()
			if hmmsearchErr != nil {
				log.Fatal(hmmsearchErr)
			} else {
				hmmsearchOutFile := outDir + "/hmmerResults.txt"
				hmmsearchOutErr := ioutil.WriteFile(hmmsearchOutFile, hmmsearchOut, 0644)
				if hmmsearchOutErr != nil {
					log.Fatal(hmmsearchOutErr)
				} else {
					hmmsearchOutStr := string(hmmsearchOut)
					hmmsearchRows := strings.Split(hmmsearchOutStr, "\n")[18:]
					// p(hmmsearchRows)
					for _, row := range hmmsearchRows {
						columns := strings.Fields(row)
						if len(columns) == 9 {
							hmmerTargets = append(hmmerTargets, columns[8])
						} else {
							break
						}
						// p(len(strings.Fields(row)), row)
					}

				}
			}

		}
	}

	// }
	return
}

func combineBlastpHmmer(blastp, hmmer []string) map[string]bool {
	combinedOut := make(map[string]bool)
	bothLists := append(blastp, hmmer...)
	for _, h := range bothLists {
		// p("header", h)
		combinedOut[h] = true
	}
	// p("BLASTP:", len(blastp))
	// p("HMMER:", len(hmmer))
	return combinedOut
}

func file2SeqMap(f string) map[string]string {
	seqOut := make(map[string]string)
	if fileExists(f) {
		fastaFh, err := os.Open(f)
		if err != nil {
			log.Fatal(err)
		}
		defer fastaFh.Close()

		for record := range parseFasta(fastaFh) {
			seqOut[record.id] = record.seq
		}

	}
	return seqOut
}

func header2seq(bothMap map[string]bool, seqFile string) map[string]string {
	seqOut := make(map[string]string)
	fastaFh, err := os.Open(seqFile)
	if err != nil {
		log.Fatal(err)
	}
	defer fastaFh.Close()

	for record := range parseFasta(fastaFh) {
		if _, ok := bothMap[record.id]; ok {
			seqOut[record.id] = record.seq
			//do something here
		}

	}
	return seqOut
}

func writeSeqMap(seqIn map[string]string, outDir, outName string) string {
	outFile := outDir + "/" + outName + ".fa"
	out, outErr := os.Create(outFile)
	if outErr != nil {
		log.Fatal(outErr)
	}
	for header, sequence := range seqIn {
		_, headerErr := out.WriteString(">" + header + "\n")
		if headerErr != nil {
			out.Close()
			log.Fatal(headerErr)
		} else {
			_, seqErr := out.WriteString(sequence + "\n")
			if seqErr != nil {
				out.Close()
				log.Fatal(seqErr)
			}
		}
	}
	closeErr := out.Close()
	if closeErr != nil {
		log.Fatal(closeErr)
	}
	return outFile
}

func concatFasta(f1, f2, outName string) (outFile string) {
	if fileExists(f1) && fileExists(f2) {
		seqMap1 := file2SeqMap(f1)
		seqMap2 := file2SeqMap(f2)
		for h, s := range seqMap1 {
			seqMap2[h] = s
		}
		outFile = writeSeqMap(seqMap2, outDir, outName)

	}
	return
}

func signalP(pepSeq map[string]string, signalPath, outDir string) (signalpOutStr string) {

	pepSeqFile := writeSeqMap(pepSeq, outDir, "preSignalP")
	signalPstr := signalPath + " -gff3 -prefix signalPout -fasta " + pepSeqFile
	signalPCmd := exec.Command("sh", "-c", signalPstr)
	timeStatus(signalPstr)
	// p("running", signalPstr)
	signalPOut, signalPErr := signalPCmd.Output()
	if signalPErr != nil {
		log.Fatal(signalPErr)
	} else {
		signalpOutStr = string(signalPOut)
	}

	return
}
func cleanHeader(originalMap map[string]string) map[string]string {
	cleanHeaderMap := make(map[string]string)
	reg, err := regexp.Compile("[^A-Za-z0-9 _  : -]+")
	if err != nil {
		log.Fatal(err)
	}

	for header := range originalMap {

		underScoreHeader := reg.ReplaceAllString(header, "_")
		// p("header", header)
		// p("underScoreHeader", underScoreHeader)
		cleanHeaderMap[underScoreHeader] = header
	}
	return cleanHeaderMap

}

// the following function should return two structures,
// 1. Header and sequence
// 2. Header and and CysPattern

func subsetSeqMap(subset map[string]int, originalMap map[string]string) (subMap, patternMap map[string]string) {

	subMap = make(map[string]string)
	patternMap = make(map[string]string)
	p("subset:", len(subset))
	p("originalMap:", len(originalMap))
	ogHeaderMap := cleanHeader(originalMap)
	p("ogHeaderMap:", len(ogHeaderMap))
	for header, cleavageSite := range subset {
		ogHeader := ogHeaderMap[header]
		if _, ok := originalMap[ogHeader]; ok {
			ogSeq := originalMap[ogHeader]
			ogMature := ogSeq[cleavageSite:]
			if len(ogSeq) < 200 && strings.Count(ogMature, "C") > 5 {
				subMap[ogHeader] = ogMature
				patternMap[ogHeader] = cysPattern(ogMature)
			}
		} else {
			p("error:", header)
		}
	}
	return
}

func signalpGffList(gffFile string) map[string]int {
	signalMap := make(map[string]int)
	signalPfileName, signalPfileError := os.Open(gffFile)
	if signalPfileError != nil {
		log.Fatal(signalPfileError)
	}
	signalPfile := csv.NewReader(signalPfileName)
	signalPfile.Comma = '\t'
	signalPfile.Comment = '#'
	withCleavage, withCleavageErr := signalPfile.ReadAll()
	if withCleavageErr != nil {
		log.Fatal(withCleavageErr)
	}

	for _, row := range withCleavage {
		header := row[0]
		cleaveSite, cleaveErr := strconv.Atoi(row[4])
		if cleaveErr != nil {
			log.Fatal(cleaveErr)
		} else {
			signalMap[header] = cleaveSite
		}

	}
	return signalMap
}

func cysPattern(aaSeq string) (pattern string) {
	var currentCysPos int
	var firstCysFound bool
	for i := 0; i < len(aaSeq); i++ {
		aa := string(aaSeq[i])
		if aa == "C" {

			if firstCysFound {
				if i-currentCysPos == 1 {
					pattern += "C"
				} else if i-currentCysPos == 2 {
					pattern += "XC"
				} else {
					pattern += "-C"
				}

			} else {
				firstCysFound = true
				pattern += "C"
			}
			currentCysPos = i
		} // end checkCys

	} // end for

	return
}

func silix(pepFileName string) (silixResults string) {

	blastp(pepFileName, pepFileName, outDir, "blastallResults", false, 1e-3)
	if willRun("silix") {
		silixStr := "silix " + pepFileName + " " + outDir + "/blastallResults.txt"
		silixCmd := exec.Command("sh", "-c", silixStr)
		timeStatus(silixStr)
		// p("running", silixStr)
		silixOut, silixErr := silixCmd.Output()
		if silixErr != nil {
			log.Fatal(silixErr)
		} else {
			silixResults = string(silixOut)
		}

	}
	return
}

func runAll(inPep string, n int) (outPep string, signalPseq, maturePattern map[string]string) {

	// var numPep int
	var bothResults map[string]bool

	signalPseq = make(map[string]string)
	maturePattern = make(map[string]string)

	// var firstRoundFinalist string
	blastResults := blastp(queryPep, inPep, outDir, "blastpResults", true, 1e-3)

	hmmerResults := hmmer(queryPep, inPep, outDir, n)

	bothResults = combineBlastpHmmer(blastResults, hmmerResults)
	blastHmmerSeqs := header2seq(bothResults, queryPep) // maybe merge combineBlastpHmmer to go here and return map[string][string]
	signalPgff := "signalPout.gff3"                     // should potentiall move the signalP gffFile
	if signalPath != "" {
		var signalpResults string

		if fileExists(signalPgff) == false {
			p(signalPgff, "doesn't exist, will generate it now")
			signalpResults = signalP(blastHmmerSeqs, signalPath, outDir)
			p(signalpResults)
			newGFF := outDir + "/signalPout_" + strconv.Itoa(n) + ".gff3"
			e := os.Rename(signalPgff, newGFF)
			if e != nil {
				log.Fatal(e)
			}
			signalPheaders := signalpGffList(newGFF)
			p(len(signalPheaders), len(blastHmmerSeqs))
			signalPseq, maturePattern = subsetSeqMap(signalPheaders, blastHmmerSeqs)
			outPep = writeSeqMap(signalPseq, outDir, "round_"+strconv.Itoa(n))
			// numPep = len(signalPseq)
			//  Here's a thought, instead of returning numPep, return map with header and pattern

		}

	} else if fileExists(signalPgff) {

		signalPheaders := signalpGffList(signalPgff)
		p(len(signalPheaders), len(blastHmmerSeqs))
		signalPseq, maturePattern = subsetSeqMap(signalPheaders, blastHmmerSeqs)
		outPep = writeSeqMap(signalPseq, outDir, "round_"+strconv.Itoa(n))
		// numPep = len(signalPseq)
	} else {
		p("Skipping signalP")
		writeSeqMap(blastHmmerSeqs, outDir, "noSignalP"+strconv.Itoa(n))
	}

	// p(time.Now().Format("15:04:05 Mon Jan-02-2006"))
	return outPep, signalPseq, maturePattern
}
func main() {
	ickTally := [2]int{}
	var numICK int
	var currentFinalist string
	var pepMap, patternMap map[string]string
	var silixData [][]string

	// pepMap = make(map[string]string)
	// patternMap = make(map[string]string)

	// startingIter := 0
	// currentOutName := "round_" + strconv.Itoa(startingIter)

	dependencies := []string{"blastp", "makeblastdb", "hmmbuild", "hmmsearch", "mafft"}
	if checkDependencies(dependencies) {
		if signalPath == "" && iterate > 1 {
			p("Please provide the path for signalp using -signalp to use -n for >1 iterations")
			p("Setting n iterations to 1")
			iterate = 1
		}

		for i := 0; i < iterate; i++ {
			if i == 0 {
				currentFinalist, pepMap, patternMap = runAll(ickDB, i)
				timeStatus("writing to: " + currentFinalist)
			} else {
				newDB := concatFasta(ickDB, currentFinalist, "newDB_"+strconv.Itoa(i))
				currentFinalist, pepMap, patternMap = runAll(newDB, i)
				numICK = len(pepMap)
				if numICK > ickTally[0] {
					ickTally[0] = numICK
					ickTally[1] = 0
				} else {
					ickTally[1]++
				}
				if ickTally[1] > 4 {
					break
				}

			}

		}
		if runSilix {

			p(patternMap)
			timeStatus(currentFinalist)
			p(fileExists(currentFinalist))
			if currentFinalist != "" && fileExists(currentFinalist) {
				silixResults := silix(currentFinalist)
				silixOutFile := outDir + "/cluck.csv"
				silixRows := splitLines(silixResults)
				for _, line := range silixRows {
					row := strings.Fields(line)
					header := row[1]
					if pattern, ok := patternMap[header]; ok {
						row = append(row, pattern)
						if maturePep, ok := pepMap[header]; ok {
							row = append(row, maturePep)
						}
						silixData = append(silixData, row)
						//do something here
					} else {
						p(header, "not found")
					}

				}

				silixOutErr := ioutil.WriteFile(silixOutFile, []byte(silixResults), 0644)
				if silixOutErr != nil {
					log.Fatal(silixOutErr)
				}

				silixCSV, silixCSVerr := os.Create(silixOutFile)
				logFatalErr("Cannot create file", silixCSVerr)
				defer silixCSV.Close()

				writer := csv.NewWriter(silixCSV)
				defer writer.Flush()

				for _, value := range silixData {
					writeErr := writer.Write(value)
					logFatalErr("Cannot write to file", writeErr)
				}
				// p(silixResults)
			}
		}
	}

}
