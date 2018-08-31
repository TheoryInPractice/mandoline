
SMALL = P2 P3 P4 P5 C4 C5 K4 K5 K3,3 S3 S4
LARGE = S5 S6 K4,4 K5,5 P6 P7 S3,2
SMALL_FILES = $(SMALL:%=example-graphs/%.dag)
LARGE_FILES = $(LARGE:%=example-graphs/%.dag)

small: ${SMALL_FILES}
large: ${LARGE_FILES}

all: small

%.dag: %.txt
	./build_counting_dag.py --quiet $< -o $@

clean:
	rm example-graphs/*.dag
