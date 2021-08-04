CC          = gcc
CFLAGS      = -Wall -O2
LDFLAGS     = -lz
prefix      = /usr/local
exec_prefix = $(prefix)/bin

src = $(wildcard *.c)
obj = $(src:.c=.o)

udini: $(obj)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) udini

.PHONY: install
install:
	cp udini $(exec_prefix)







