package tools

import (
	"container/list"
	"image"
	"image/color"
	"image/draw"
	"math"
)

type point struct {
	X, Y float64
}

type imgPointer struct {
	*image.RGBA
}

type points []point

type segment struct {
	old, new point
}

type segments []segment

type character struct {
	strokes   segments
	thickness int
	c         color.RGBA
}

type characters []character

type cBox struct {
	boundry, interior, filled color.RGBA
}

type circle struct {
	radius    float64
	center    point
	thickness int
}

const (
	grid_size            = 100000
	par_space_upperbound = 2.0 * math.Pi
	par_space_lowerbound = float64(0)
	factor_X             = 20
	factor_Y             = 20
	lowerLeft_X          = 0
	lowerLeft_Y          = 0
	upperRight_X         = 1024
	upperRight_Y         = 1024
	theta                = (-1.0) * math.Pi
)

var (
	factor, displacement  point      = point{factor_X, factor_Y}, point{upperRight_X / 2, upperRight_Y * 0.4}
	lowerLeft, upperRight point      = point{lowerLeft_X, lowerLeft_Y}, point{upperRight_X, upperRight_Y}
	par_space             point      = point{par_space_lowerbound, par_space_upperbound}
	red, white, blue      color.RGBA = color.RGBA{255, 0, 0, 255}, color.RGBA{255, 255, 255, 255}, color.RGBA{0, 0, 255, 255}
	T                     character  = character{segments{segment{point{300, 800}, point{10000, 800}},
		segment{point{450, 905}, point{450, 950}}}, 1, blue}
	unitEast, unitWest, unitSouth, unitNorth point = point{1, 0}, point{-1, 0}, point{0, -1}, point{0, 1}
)

func drr() imgPointer {
	grid, target := makeData(grid_size, par_space)
	target.maps(grid)
	return newCanvas(lowerLeft, upperRight, white).
		markAll(target).
		fill(point{512, 512}, cBox{red, white, red}).
		hLines(190, 250, 800, 1, blue).
		vLines(800, 900, 220, 1, blue).
		hLines(190, 250, 900, 1, blue).
		lines(point{310, 900}, point{330, 800}, 1, blue).
		lines(point{330, 800}, point{350, 900}, 1, blue).
		lines(point{350, 900}, point{370, 800}, 1, blue).
		lines(point{370, 800}, point{390, 900}, 1, blue).
		vLines(830, 835, 405, 1, blue).
		vLines(855, 900, 405, 1, blue).
		drawArch(circle{12, point{433, 862}, 1}, math.Pi/float64(2), math.Pi*float64(3)/float64(2)+0.05, blue).
		drawArch(circle{12, point{433, 886}, 1}, math.Pi*float64(3)/float64(2), math.Pi/float64(2)+0.2, blue).
		drawArch(circle{12, point{463, 862}, 1}, math.Pi/float64(2), math.Pi*float64(3)/float64(2)+0.05, blue).
		drawArch(circle{12, point{463, 886}, 1}, math.Pi*float64(3)/float64(2), math.Pi/float64(2)+0.1, blue).
		lines(point{510, 800}, point{540, 850}, 1, blue).
		lines(point{540, 850}, point{570, 800}, 1, blue).
		lines(point{540, 850}, point{540, 900}, 1, blue).
		drawCircle(circle{24, point{575, 875}, 1}, blue).
		vLines(850, 885, 610, 1, blue).
		vLines(850, 900, 640, 1, blue).
		drawArch(circle{15, point{625, 885}, 1}, 0, math.Pi, blue).
		vLines(800, 900, 700, 1, blue).
		hLines(700, 760, 900, 1, blue).
		vLines(899, 902, 770, 2, blue).
		vLines(800, 900, 800, 1, blue).
		hLines(800, 860, 900, 1, blue).
		vLines(899, 902, 870, 2, blue)

}

func makeData(grid_size int, par_space point) ([]float64, points) {
	grid := make([]float64, grid_size)
	for i := 0; i < grid_size; i++ {
		grid[i] = par_space.X + float64(i)/float64(grid_size)*(par_space.Y-par_space.X)
	}
	grid = append(grid, par_space_upperbound)
	var target points
	//fmt.Println(grid)
	//fmt.Println(target)
	return grid, target
}

func (p *point) stretch(factor, displacement point) *point {
	p.X = p.X*factor.X + displacement.X
	p.Y = p.Y*factor.Y + displacement.Y
	//fmt.Println(p)
	return p
}

func (p *point) calc(t float64) *point {
	p.X = 16 * math.Sin(t) * math.Sin(t) * math.Sin(t)
	p.Y = 13*math.Cos(t) - 5*math.Cos(2*t) - 2*math.Cos(3*t) - math.Cos(4*t)
	//fmt.Println(p)
	return p
}

func (p *point) rotate(theta float64) *point {
	p.X = p.X*math.Cos(theta) - p.Y*math.Sin(theta)
	p.Y = p.X*math.Sin(theta) + p.Y*math.Cos(theta)
	return p
}

func (target *points) maps(grid []float64) {
	for _, t := range grid {
		var p point
		p.calc(t).rotate(theta).stretch(factor, displacement)
		*target = append(*target, p)
	}
}

func (img imgPointer) mark(p point, c color.RGBA) {
	img.SetRGBA(int(p.X), int(p.Y), c)
}

func (heart imgPointer) markAll(target []point) imgPointer {
	for _, p := range target {
		heart.mark(p, red)
	}
	return heart
}

func (heart imgPointer) drawCircle(circle circle, c color.RGBA) imgPointer {
	for x := int(circle.center.X) - int(circle.radius) - 2*circle.thickness; x < int(circle.center.X)+int(circle.radius)+2*circle.thickness; x++ {
		for y := int(circle.center.Y) - int(circle.radius) - 2*circle.thickness; y < int(circle.center.Y)+int(circle.radius)+2*circle.thickness; y++ {
			p := point{float64(x), float64(y)}
			if p.dist(circle.center) <= circle.radius+float64(circle.thickness) && p.dist(circle.center) >= circle.radius-float64(circle.thickness) {
				heart.SetRGBA(x, y, c)
			}
		}
	}
	return heart
}

func (heart imgPointer) drawArch(circle circle, start, end float64, c color.RGBA) imgPointer {
	for x := int(circle.center.X) - int(circle.radius) - 2*circle.thickness; x < int(circle.center.X)+int(circle.radius)+2*circle.thickness; x++ {
		for y := int(circle.center.Y) - int(circle.radius) - 2*circle.thickness; y < int(circle.center.Y)+int(circle.radius)+2*circle.thickness; y++ {
			p := point{float64(x), float64(y)}
			if p.dist(circle.center) <= circle.radius+float64(circle.thickness) && p.dist(circle.center) >= circle.radius-float64(circle.thickness) && p.isWithin(circle, start, end) {
				heart.SetRGBA(x, y, c)
			}
		}
	}
	return heart
}

func newCanvas(lowerLeft, upperRight point, c color.RGBA) imgPointer {
	img := image.NewRGBA(image.Rect(
		int(lowerLeft.X), int(lowerLeft.Y),
		int(upperRight_X), int(upperRight_Y)))
	draw.Draw(img, img.Bounds(), &image.Uniform{c}, image.Point{0, 0}, draw.Src)
	return imgPointer{img}
}

func (heart imgPointer) hLine(x_old, x_new, y int, c color.RGBA) imgPointer {
	for ; x_old <= x_new; x_old++ {
		heart.SetRGBA(x_old, y, c)
	}
	return heart
}

func (heart imgPointer) vLine(y_old, y_new, x int, c color.RGBA) imgPointer {
	for ; y_old <= y_new; y_old++ {
		heart.SetRGBA(x, y_old, c)
	}
	return heart
}

func (heart imgPointer) hLines(x_old, x_new, y, thickness int, c color.RGBA) imgPointer {
	for i := 0; i < 1+2*thickness; i++ {
		heart.hLine(x_old, x_new, y-thickness+i, c)
	}
	return heart
}

func (heart imgPointer) vLines(y_old, y_new, x, thickness int, c color.RGBA) imgPointer {
	for i := 0; i < 1+2*thickness; i++ {
		heart.vLine(y_old, y_new, x-thickness+i, c)
	}
	return heart
}

func (heart imgPointer) concatenateHorizontalLines(old, new point, thickness int, c color.RGBA) imgPointer {
	var cutoffs points = calcCutOffs(old, new, int(math.Abs(new.Y-old.Y)))
	for i := range cutoffs {
		heart.hLines(int(cutoffs[i].X), int(cutoffs[i+1].X), int(cutoffs[i].Y), thickness, c)
	}
	return heart
}

func (heart imgPointer) concatenateVerticalLines(old, new point, thickness int, c color.RGBA) imgPointer {
	var cutoffs points = calcCutOffs(old, new, int(math.Abs(new.X-old.X)))
	for i := range cutoffs {
		heart.vLines(int(cutoffs[i].Y), int(cutoffs[i+1].Y), int(cutoffs[i].X), thickness, c)
	}
	return heart
}

func (heart imgPointer) drawFlatLine(old, new point, thickness int, c color.RGBA) imgPointer {
	old, new = reorderX(old, new)
	slope := (new.Y - old.Y) / (new.X - old.X)
	for i := 0; i <= int(new.X)-int(old.X); i++ {
		heart.vLine(int(old.Y+slope*float64(i))-thickness, int(old.Y+slope*float64(i))+thickness, int(old.X)+i, c)
	}
	return heart
}

func (heart imgPointer) drawSteepLine(old, new point, thickness int, c color.RGBA) imgPointer {
	old, new = reorderY(old, new)
	slope := (new.Y - old.Y) / (new.X - old.X)
	for i := 0; i <= int(new.Y)-int(old.Y); i++ {
		heart.hLine(int(old.X)+int(float64(i)/slope)-thickness, int(old.X)+int(float64(i)/slope)+thickness, int(old.Y)+i, c)
	}
	return heart
}

func reorderY(old, new point) (point, point) {
	if old.Y >= new.Y {
		temp := old
		old = new
		new = temp
		return old, new
	} else {
		return old, new
	}
}

func reorderX(old, new point) (point, point) {
	if old.X >= new.X {
		temp := old
		old = new
		new = temp
		return old, new
	} else {
		return old, new
	}
}

func isLessThanPiOverTwo(old, new point) bool {
	return math.Abs((new.Y-old.Y)/(new.X-old.X)) >= 1.0
}

func isVertical(old, new point) bool {
	return old.X == new.X
}

func isHorizontal(old, new point) bool {
	return old.Y == new.Y
}

func judge(old, new point) string {
	vertical, horizontal, piOverTwo := isVertical(old, new), isHorizontal(old, new), isLessThanPiOverTwo(old, new)
	switch {
	case vertical:
		return "v"
	case horizontal:
		return "h"
	case piOverTwo:
		return "iv"
	default:
		return "ih"
	}
}

func calcCutOffs(old, new point, num int) points {
	cutoffs, diff := make([]point, num), new.add(old.scalarMultiply(float64(-1)))
	for i := range cutoffs {
		cutoffs[i] = old.add(diff.scalarMultiply(float64(i+1) / float64(num+1)))
	}
	return append(cutoffs, new)
}

func (u point) add(v point) point {
	return point{u.X + v.X, u.Y + v.Y}
}

func (u point) scalarMultiply(c float64) point {
	return point{c * u.X, c * u.Y}
}

func (u point) equals(v point) bool {
	return (u.X == v.X) && (u.Y == v.Y)
}

func (u point) dist(v point) float64 {
	return math.Sqrt(math.Pow(u.X-v.X, 2) + math.Pow(u.Y-v.Y, 2))
}

func (u point) reflect() point {
	x, y := u.Y, u.X
	return point{x, y}
}

func (p point) isWithin(circle circle, start, end float64) bool {
	pStart, pEnd := point{circle.center.X + circle.radius*math.Cos(start), circle.center.Y + circle.radius*math.Sin(start)},
		point{circle.center.X + circle.radius*math.Cos(end), circle.center.Y + circle.radius*math.Sin(end)}
	pPos, pStartPos, pEndPos := p.quadrant(circle.center), pStart.quadrant(circle.center), pEnd.quadrant(circle.center)
	if start > end {
		return (pPos >= pStartPos) || (pPos <= pEndPos)
	} else {
		if pPos == pStartPos {
			pStartPos = testStartAngle(circle, p, start, pPos)
		}
		if pPos == pEndPos {
			pEndPos = testEndAngle(circle, p, end, pPos)
		}
		return (pPos >= pStartPos) && (pPos <= pEndPos)
	}

}

func testStartAngle(circle circle, p point, start, pPos float64) float64 {
	sineP, sineStart := (p.Y-circle.center.Y)/p.dist(circle.center), math.Sin(start)
	switch pPos {
	case 2, 8:
		if sineP >= sineStart {
			return pPos
		} else {
			return pPos + 0.5
		}
	case 4, 6:
		if sineP <= sineStart {
			return pPos
		} else {
			return pPos + 0.5
		}
	default:
		return pPos
	}
}

func testEndAngle(circle circle, p point, end, pPos float64) float64 {
	sineP, sineEnd := (p.Y-circle.center.Y)/p.dist(circle.center), math.Sin(end)
	switch pPos {
	case 2, 8:
		if sineP <= sineEnd {
			return pPos
		} else {
			return pPos + 0.5
		}
	case 4, 6:
		if sineP >= sineEnd {
			return pPos
		} else {
			return pPos + 0.5
		}
	default:
		return pPos
	}
}

func (p point) quadrant(origin point) float64 {
	xPositive, yPositive := p.X > origin.X, p.Y > origin.Y
	xNegative, yNegative := p.X < origin.X, p.Y < origin.Y
	yAxis, xAxis := p.X == origin.X, p.Y == origin.Y
	switch {
	case xPositive && yPositive:
		return 2
	case xPositive && yNegative:
		return 8
	case xNegative && yNegative:
		return 6
	case xNegative && yPositive:
		return 4
	case yAxis && yPositive:
		return 3
	case yAxis && yNegative:
		return 7
	case xAxis && xPositive:
		return 1
	case xAxis && xNegative:
		return 5
	default:
		return 0
	}
}

func (heart imgPointer) lines(old, new point, thickness int, c color.RGBA) imgPointer {
	switch result := judge(old, new); result {
	case "v":
		heart.vLines(int(old.Y), int(new.Y), int(old.X), thickness, c)
	case "h":
		heart.hLines(int(old.X), int(new.X), int(old.Y), thickness, c)
	case "iv":
		heart.drawSteepLine(old, new, thickness, c)
	default:
		heart.drawFlatLine(old, new, thickness, c)
	}
	return heart
}

func (heart imgPointer) fillLine(seed point, box cBox) imgPointer {
	s := findSegment(seed, box, heart)
	heart.fillSegment(s, box.filled)
	return heart
}

func (heart imgPointer) naiveFillConvex(seed point, box cBox) imgPointer {
	heart.fillLine(seed, box)
	newSeedUp, newSeedDown := seed.add(point{0, 1}), seed.add(point{0, -1})
	for heart.notFilled(newSeedUp, box.interior) || heart.notFilled(newSeedDown, box.interior) {
		if heart.notFilled(newSeedUp, box.interior) {
			heart.fillLine(newSeedUp, box)
			newSeedUp = newSeedUp.add(point{0, 1})
		}
		if heart.notFilled(newSeedDown, box.interior) {
			heart.fillLine(newSeedDown, box)
			newSeedDown = newSeedDown.add(point{0, -1})
		}
	}
	return heart
}

func (heart imgPointer) notFilled(p point, interior color.RGBA) bool {
	return heart.RGBAAt(int(p.X), int(p.Y)) == interior
}

func (heart imgPointer) fill(seed point, box cBox) imgPointer {
	frontier := list.New()
	frontier.PushBack(seed)
	heart.SetRGBA(int(seed.X), int(seed.Y), box.filled)
	for frontier.Len() != 0 {
		filledPoint := frontier.Remove(frontier.Front()).(point)
		east, west, south, north := filledPoint.add(unitEast), filledPoint.add(unitWest), filledPoint.add(unitSouth), filledPoint.add(unitNorth)
		if heart.notFilled(east, box.interior) {
			frontier.PushBack(east)
			heart.SetRGBA(int(east.X), int(east.Y), box.filled)
		}
		if heart.notFilled(west, box.interior) {
			frontier.PushBack(west)
			heart.SetRGBA(int(west.X), int(west.Y), box.filled)
		}
		if heart.notFilled(south, box.interior) {
			frontier.PushBack(south)
			heart.SetRGBA(int(south.X), int(south.Y), box.filled)
		}
		if heart.notFilled(north, box.interior) {
			frontier.PushBack(north)
			heart.SetRGBA(int(north.X), int(north.Y), box.filled)
		}
	}
	return heart
}

func (heart imgPointer) fillSegment(s segment, c color.RGBA) imgPointer {
	heart.hLine(int(s.old.X), int(s.new.X), int(s.old.Y), c)
	return heart
}

func findSegment(seed point, box cBox, heart imgPointer) segment {
	var old, new point
	x, y := int(seed.X), int(seed.Y)
	for heart.RGBAAt(x, y) == box.interior {
		x -= 1
	}
	old = point{float64(x + 1), float64(y)}
	x = int(seed.X)
	for heart.RGBAAt(x, y) == box.interior {
		x += 1
	}
	new = point{float64(x - 1), float64(y)}
	return segment{old, new}

}

func (heart imgPointer) fillSegments(unfilledSegments segments, c color.RGBA) imgPointer {
	for _, s := range unfilledSegments {
		heart.fillSegment(s, c)
	}
	return heart
}

func findSegments(filledSegment segment, box cBox, heart imgPointer) segments {
	unfilledSegments := segments{}
	unfilledSegments.findSegmentsOneDirection(filledSegment, box, heart, 1).
		findSegmentsOneDirection(filledSegment, box, heart, -1)
	return unfilledSegments
}

func (s *segments) findSegmentsOneDirection(filledSegment segment, box cBox, heart imgPointer, direction int) *segments {
	unfilled, zeroLength := heart.findSeed(filledSegment, box, direction), filledSegment.new.equals(filledSegment.old)
	switch {
	case zeroLength && unfilled:
		*s = append(*s, getUnfilledPoint(filledSegment, direction))
		return s
	case zeroLength && (!unfilled):
		return s
	case (!zeroLength) && (!unfilled):
		left, right := split(filledSegment)
		s.findSegmentsOneDirection(left, box, heart, direction)
		s.findSegmentsOneDirection(right, box, heart, direction)
		return s
	default:
		newSegment := findSegment(getSeed(filledSegment, direction), box, heart)
		*s = append(*s, newSegment)
		exceedLeft, exceedRight := newSegment.old.X <= filledSegment.old.X, newSegment.new.X >= filledSegment.new.X
		switch {
		case exceedLeft && exceedRight:
			return s
		case exceedLeft && (!exceedRight):
			newRightFilledSegment := getRightFilledSegment(filledSegment, newSegment)
			s.findSegmentsOneDirection(newRightFilledSegment, box, heart, direction)
			return s
		case (!exceedLeft) && exceedRight:
			newLeftFilledSegment := getLeftFilledSegment(filledSegment, newSegment)
			s.findSegmentsOneDirection(newLeftFilledSegment, box, heart, direction)
			return s
		default:
			newRightFilledSegment := getRightFilledSegment(filledSegment, newSegment)
			newLeftFilledSegment := getLeftFilledSegment(filledSegment, newSegment)
			s.findSegmentsOneDirection(newRightFilledSegment, box, heart, direction)
			s.findSegmentsOneDirection(newLeftFilledSegment, box, heart, direction)
			return s
		}

	}

}

func getUnfilledPoint(s segment, direction int) segment {
	return segment{s.old.add(point{0, float64(direction)}), s.new.add(point{0, float64(direction)})}
}
func split(s segment) (left, right segment) {
	mid := getMid(s)
	return segment{s.old, mid.add(point{-1, 0})}, segment{mid, s.new}
}
func getMid(s segment) point {
	return s.new.scalarMultiply(0.5).add(s.old.scalarMultiply(0.5))
}
func getRightFilledSegment(filledSegment, newSegment segment) segment {
	return segment{newSegment.new.add(point{float64(1), 0}), filledSegment.new}
}
func getLeftFilledSegment(filledSegment, newSegment segment) segment {
	return segment{filledSegment.old, newSegment.old.add(point{float64(-1), 0})}
}

func (heart imgPointer) findSeed(filledSegment segment, box cBox, direction int) bool {
	mid := getMid(filledSegment)
	switch c := heart.RGBAAt(int(mid.X), int(mid.Y)+direction); c {
	case box.boundry:
		return false
	case box.filled:
		return false
	default:
		return true
	}
}

func getSeed(filledSegment segment, direction int) point {
	mid := getMid(filledSegment)
	return mid.add(point{0, float64(direction)})
}
