TARBALLS_DIR := tarballs
REPOS := $(wildcard packages/* tools/*)
TARBALLS := $(patsubst %,${TARBALLS_DIR}/%.tar.gz,${REPOS})

.PHONY: all clean

all: ${TARBALLS}

clean:
	rm -rf tarballs

${TARBALLS_DIR}/%.tar.gz: %
	@mkdir -p $(@D)
	tar -L --directory $< --create --gzip --file $@ .
