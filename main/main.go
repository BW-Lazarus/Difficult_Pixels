package main

import (
	"image/png"
	"os"
)

func main() {

	img := tools.drr()
	output, _ := os.Create("your_pic.png")
	png.Encode(output, img)
	output.Close()
}
