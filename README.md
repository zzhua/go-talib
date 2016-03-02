# go-talib

[![GoDoc](http://godoc.org/github.com/markcheno/go-talib?status.svg)](http://godoc.org/github.com/markcheno/go-talib) 

A pure [Go](http://golang.org/) port of [TA-Lib](http://ta-lib.org)

## Install

Install the package with:

```bash
go get github.com/markcheno/go-talib
```

Import it with:

```go
import "github.com/markcheno/go-talib"
```

and use `talib` as the package name inside the code.

## Example

```go
package main

import (
  "fmt"
  "github.com/markcheno/go-talib"
)

func main() {

  var data = []float64{201.28, 197.64, 195.78, 198.22, 201.74, 200.12, 198.55, 197.99, 196.80, 195.00,
                       197.55, 197.97, 198.97, 201.93, 200.83, 201.30, 198.64, 196.09, 197.91, 195.42,
                       197.84, 200.70, 199.93, 201.95, 201.39, 200.49, 202.63, 202.75, 204.70, 205.54}

  sma := talib.Sma(data, 5)
  
  fmt.Println(sma)
}
```

## License

MIT License  - see LICENSE for more details
