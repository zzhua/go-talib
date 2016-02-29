// Package talib is a pure Go port of <a href=http://ta-lib.org>TaLib</a> Technical Analysis Library
// Licensed under terms of MIT license (see LICENSE)
// Copyright (c) 2016 Mark Chenoweth, mark@cheno.net
package talib

import (
	"math"
)

/* Overlap Studies
TODO:
  BBANDS - Bollinger Bands
    upperband, middleband, lowerband = BBANDS(close, timeperiod=5, nbdevup=2, nbdevdn=2, matype=0)
  DEMA - Double Exponential Moving Average
    real = DEMA(close, timeperiod=30)
  HT_TRENDLINE - Hilbert Transform - Instantaneous Trendline
    real = HT_TRENDLINE(close)
  KAMA - Kaufman Adaptive Moving Average
    real = KAMA(close, timeperiod=30)
  MA - Moving average
    real = MA(close, timeperiod=30, matype=0)
  MAMA - MESA Adaptive Moving Average
    mama, fama = MAMA(close, fastlimit=0, slowlimit=0)
  MAVP - Moving average with variable period
    real = MAVP(close, periods, minperiod=2, maxperiod=30, matype=0)
  MIDPOINT - MidPoint over period
    real = MIDPOINT(close, timeperiod=14)
  MIDPRICE - Midpoint Price over period
    real = MIDPRICE(high, low, timeperiod=14)
  SAR - Parabolic SAR
    real = SAR(high, low, acceleration=0, maximum=0)
  SAREXT - Parabolic SAR - Extended
    real = SAREXT(high, low, startvalue=0, offsetonreverse=0, accelerationinitlong=0, accelerationlong=0, accelerationmaxlong=0, accelerationinitshort=0, accelerationshort=0, accelerationmaxshort=0)
  T3 - Triple Exponential Moving Average (T3)
    real = T3(close, timeperiod=5, vfactor=0)
  TEMA - Triple Exponential Moving Average
    real = TEMA(close, timeperiod=30)
  TRIMA - Triangular Moving Average
    real = TRIMA(close, timeperiod=30)
  WMA - Weighted Moving Average
    real = WMA(close, timeperiod=30)
*/

// Ema - Exponential Moving Average
func Ema(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))
	optInK1 := (2.0 / (float64(optInTimePeriod) + 1))

	lookbackTotal := optInTimePeriod - 1
	startIdx := lookbackTotal
	today := startIdx - lookbackTotal
	i := optInTimePeriod
	tempReal := 0.0
	for i > 0 {
		tempReal += inReal[today]
		today++
		i--
	}
	prevMA := tempReal / float64(optInTimePeriod)
	for today <= startIdx {
		prevMA = ((inReal[today] - prevMA) * optInK1) + prevMA
		today++
	}
	outReal[startIdx] = prevMA
	outIdx := startIdx + 1

	for today < len(outReal) {
		prevMA = ((inReal[today] - prevMA) * optInK1) + prevMA
		outReal[outIdx] = prevMA
		today++
		outIdx++
	}

	return outReal
}

// Sma - Simple Moving Average
func Sma(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	lookbackTotal := optInTimePeriod - 1
	startIdx := lookbackTotal
	periodTotal := 0.0
	trailingIdx := startIdx - lookbackTotal
	i := trailingIdx
	if optInTimePeriod > 1 {
		for i < startIdx {
			periodTotal += inReal[i]
			i++
		}
	}
	outIdx := startIdx
	for ok := true; ok; {
		periodTotal += inReal[i]
		tempReal := periodTotal
		periodTotal -= inReal[trailingIdx]
		outReal[outIdx] = tempReal / float64(optInTimePeriod)
		trailingIdx++
		i++
		outIdx++
		ok = i < len(outReal)
	}

	return outReal
}

/* Momentum Indicators
TODO:
  APO - Absolute Price Oscillator
    real = APO(close, fastperiod=12, slowperiod=26, matype=0)
  AROON - Aroon
    aroondown, aroonup = AROON(high, low, timeperiod=14)
  AROONOSC - Aroon Oscillator
    real = AROONOSC(high, low, timeperiod=14)
  BOP - Balance Of Power
    real = BOP(open, high, low, close)
  CCI - Commodity Channel Index
    real = CCI(high, low, close, timeperiod=14)
  CMO - Chande Momentum Oscillator
    real = CMO(close, timeperiod=14)
  DX - Directional Movement Index
    real = DX(high, low, close, timeperiod=14)
  MACD - Moving Average Convergence/Divergence
    macd, macdsignal, macdhist = MACD(close, fastperiod=12, slowperiod=26, signalperiod=9)
  MACDEXT - MACD with controllable MA type
    macd, macdsignal, macdhist = MACDEXT(close, fastperiod=12, fastmatype=0, slowperiod=26, slowmatype=0, signalperiod=9, signalmatype=0)
  MACDFIX - Moving Average Convergence/Divergence Fix 12/26
    macd, macdsignal, macdhist = MACDFIX(close, signalperiod=9)
  MFI - Money Flow Index
    real = MFI(high, low, close, volume, timeperiod=14)
  MINUS_DI - Minus Directional Indicator
    real = MINUS_DI(high, low, close, timeperiod=14)
  MINUS_DM - Minus Directional Movement
    real = MINUS_DM(high, low, timeperiod=14)
  MOM - Momentum
    real = MOM(close, timeperiod=10)
  PLUS_DI - Plus Directional Indicator
    real = PLUS_DI(high, low, close, timeperiod=14)
  PLUS_DM - Plus Directional Movement
    real = PLUS_DM(high, low, timeperiod=14)
  PPO - Percentage Price Oscillator
    real = PPO(close, fastperiod=12, slowperiod=26, matype=0)
  ROC - Rate of change : ((price/prevPrice)-1)*100
    real = ROC(close, timeperiod=10)
  ROCR - Rate of change ratio: (price/prevPrice)
    real = ROCR(close, timeperiod=10)
  ROCR100 - Rate of change ratio 100 scale: (price/prevPrice)*100
    real = ROCR100(close, timeperiod=10)
  STOCH - Stochastic
    slowk, slowd = STOCH(high, low, close, fastk_period=5, slowk_period=3, slowk_matype=0, slowd_period=3, slowd_matype=0)
  STOCHF - Stochastic Fast
    fastk, fastd = STOCHF(high, low, close, fastk_period=5, fastd_period=3, fastd_matype=0)
  STOCHRSI - Stochastic Relative Strength Index
    fastk, fastd = STOCHRSI(close, timeperiod=14, fastk_period=5, fastd_period=3, fastd_matype=0)
  TRIX - 1-day Rate-Of-Change (ROC) of a Triple Smooth EMA
    real = TRIX(close, timeperiod=30)
  ULTOSC - Ultimate Oscillator
    real = ULTOSC(high, low, close, timeperiod1=7, timeperiod2=14, timeperiod3=28)
*/

// Adx - Average Directional Movement Index
func Adx(inHigh []float64, inLow []float64, inClose []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inClose))

	optInTimePeriodF := float64(optInTimePeriod)
	lookbackTotal := (2 * optInTimePeriod) - 1
	startIdx := lookbackTotal
	outIdx := optInTimePeriod
	prevMinusDM := 0.0
	prevPlusDM := 0.0
	prevTR := 0.0
	today := startIdx - lookbackTotal
	prevHigh := inHigh[today]
	prevLow := inLow[today]
	prevClose := inClose[today]
	for i := optInTimePeriod - 1; i > 0; i-- {
		today++
		tempReal := inHigh[today]
		diffP := tempReal - prevHigh
		prevHigh = tempReal
		tempReal = inLow[today]
		diffM := prevLow - tempReal
		prevLow = tempReal
		if (diffM > 0) && (diffP < diffM) {
			prevMinusDM += diffM
		} else if (diffP > 0) && (diffP > diffM) {
			prevPlusDM += diffP
		}
		tempReal = prevHigh - prevLow
		tempReal2 := math.Abs(prevHigh - prevClose)
		if tempReal2 > tempReal {
			tempReal = tempReal2
		}
		tempReal2 = math.Abs(prevLow - prevClose)
		if tempReal2 > tempReal {
			tempReal = tempReal2
		}

		prevTR += tempReal
		prevClose = inClose[today]
	}
	sumDX := 0.0
	for i := optInTimePeriod; i > 0; i-- {
		today++
		tempReal := inHigh[today]
		diffP := tempReal - prevHigh
		prevHigh = tempReal
		tempReal = inLow[today]
		diffM := prevLow - tempReal
		prevLow = tempReal
		prevMinusDM -= prevMinusDM / optInTimePeriodF
		prevPlusDM -= prevPlusDM / optInTimePeriodF
		if (diffM > 0) && (diffP < diffM) {
			prevMinusDM += diffM
		} else if (diffP > 0) && (diffP > diffM) {
			prevPlusDM += diffP
		}
		tempReal = prevHigh - prevLow
		tempReal2 := math.Abs(prevHigh - prevClose)
		if tempReal2 > tempReal {
			tempReal = tempReal2
		}
		tempReal2 = math.Abs(prevLow - prevClose)
		if tempReal2 > tempReal {
			tempReal = tempReal2
		}

		prevTR = prevTR - (prevTR / optInTimePeriodF) + tempReal
		prevClose = inClose[today]
		if !(((-(0.00000000000001)) < prevTR) && (prevTR < (0.00000000000001))) {
			minusDI := (100.0 * (prevMinusDM / prevTR))
			plusDI := (100.0 * (prevPlusDM / prevTR))
			tempReal = minusDI + plusDI
			if !(((-(0.00000000000001)) < tempReal) && (tempReal < (0.00000000000001))) {
				sumDX += (100.0 * (math.Abs(minusDI-plusDI) / tempReal))
			}
		}
	}
	prevADX := (sumDX / optInTimePeriodF)

	outReal[startIdx] = prevADX
	outIdx = startIdx + 1
	today++
	for today < len(inClose) {
		tempReal := inHigh[today]
		diffP := tempReal - prevHigh
		prevHigh = tempReal
		tempReal = inLow[today]
		diffM := prevLow - tempReal
		prevLow = tempReal
		prevMinusDM -= prevMinusDM / optInTimePeriodF
		prevPlusDM -= prevPlusDM / optInTimePeriodF
		if (diffM > 0) && (diffP < diffM) {
			prevMinusDM += diffM
		} else if (diffP > 0) && (diffP > diffM) {
			prevPlusDM += diffP
		}
		tempReal = prevHigh - prevLow
		tempReal2 := math.Abs(prevHigh - prevClose)
		if tempReal2 > tempReal {
			tempReal = tempReal2
		}
		tempReal2 = math.Abs(prevLow - prevClose)
		if tempReal2 > tempReal {
			tempReal = tempReal2
		}

		prevTR = prevTR - (prevTR / optInTimePeriodF) + tempReal
		prevClose = inClose[today]
		if !(((-(0.00000000000001)) < prevTR) && (prevTR < (0.00000000000001))) {
			minusDI := (100.0 * (prevMinusDM / prevTR))
			plusDI := (100.0 * (prevPlusDM / prevTR))
			tempReal = minusDI + plusDI
			if !(((-(0.00000000000001)) < tempReal) && (tempReal < (0.00000000000001))) {
				tempReal = (100.0 * (math.Abs(minusDI-plusDI) / tempReal))
				prevADX = (((prevADX * (optInTimePeriodF - 1)) + tempReal) / optInTimePeriodF)
			}
		}
		outReal[outIdx] = prevADX
		outIdx++
		today++
	}
	return outReal
}

// AdxR - Average Directional Movement Index Rating
func AdxR(inHigh []float64, inLow []float64, inClose []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inClose))
	startIdx := (2 * optInTimePeriod) - 1
	tmpadx := Adx(inHigh, inLow, inClose, optInTimePeriod)
	i := startIdx
	j := startIdx + optInTimePeriod - 1
	for outIdx := startIdx + optInTimePeriod - 1; outIdx < len(inClose); outIdx, i, j = outIdx+1, i+1, j+1 {
		outReal[outIdx] = ((tmpadx[i] + tmpadx[j]) / 2.0)
	}
	return outReal
}

// Rocp - Rate of change Percentage: (price-prevPrice)/prevPrice
func Rocp(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	if optInTimePeriod < 1 {
		return outReal
	}

	startIdx := optInTimePeriod
	outIdx := startIdx
	inIdx := startIdx
	trailingIdx := startIdx - optInTimePeriod
	for inIdx < len(outReal) {
		tempReal := inReal[trailingIdx]
		if tempReal != 0.0 {
			outReal[outIdx] = (inReal[inIdx] - tempReal) / tempReal
		} else {
			outReal[outIdx] = 0.0
		}
		trailingIdx++
		outIdx++
		inIdx++
	}

	return outReal
}

// Rsi - Relative strength index
func Rsi(inReal []float64, inTimePeriod int) []float64 {

	optInTimePeriod := float64(inTimePeriod)
	outReal := make([]float64, len(inReal))

	if optInTimePeriod < 2 {
		return outReal
	}

	// variable declarations
	tempValue1 := 0.0
	tempValue2 := 0.0
	outIdx := inTimePeriod //0
	today := 0
	prevValue := inReal[today]
	prevGain := 0.0
	prevLoss := 0.0
	today++

	for i := optInTimePeriod; i > 0; i-- {
		tempValue1 = inReal[today]
		today++
		tempValue2 = tempValue1 - prevValue
		prevValue = tempValue1
		if tempValue2 < 0 {
			prevLoss -= tempValue2
		} else {
			prevGain += tempValue2
		}
	}

	prevLoss /= optInTimePeriod
	prevGain /= optInTimePeriod

	if today > 0 {

		tempValue1 = prevGain + prevLoss
		if !((-0.00000000000001 < tempValue1) && (tempValue1 < 0.00000000000001)) {
			outReal[outIdx] = 100.0 * (prevGain / tempValue1)
		} else {
			outReal[outIdx] = 0.0
		}
		outIdx++

	} else {

		for today < 0 {
			tempValue1 = inReal[today]
			tempValue2 = tempValue1 - prevValue
			prevValue = tempValue1
			prevLoss *= (optInTimePeriod - 1)
			prevGain *= (optInTimePeriod - 1)
			if tempValue2 < 0 {
				prevLoss -= tempValue2
			} else {
				prevGain += tempValue2
			}
			prevLoss /= optInTimePeriod
			prevGain /= optInTimePeriod
			today++
		}
	}

	for today < len(inReal) {

		tempValue1 = inReal[today]
		today++
		tempValue2 = tempValue1 - prevValue
		prevValue = tempValue1
		prevLoss *= (optInTimePeriod - 1)
		prevGain *= (optInTimePeriod - 1)
		if tempValue2 < 0 {
			prevLoss -= tempValue2
		} else {
			prevGain += tempValue2
		}
		prevLoss /= optInTimePeriod
		prevGain /= optInTimePeriod
		tempValue1 = prevGain + prevLoss
		if !((-0.00000000000001 < tempValue1) && (tempValue1 < 0.00000000000001)) {
			outReal[outIdx] = 100.0 * (prevGain / tempValue1)
		} else {
			outReal[outIdx] = 0.0
		}
		outIdx++
	}

	return outReal
}

// WillR - Williams' %R
func WillR(inHigh []float64, inLow []float64, inClose []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inClose))
	nbInitialElementNeeded := (optInTimePeriod - 1)
	diff := 0.0
	outIdx := optInTimePeriod - 1
	startIdx := optInTimePeriod - 1
	today := startIdx
	trailingIdx := startIdx - nbInitialElementNeeded
	highestIdx := -1
	lowestIdx := -1
	highest := 0.0
	lowest := 0.0
	i := 0
	for today < len(inClose) {
		tmp := inLow[today]
		if lowestIdx < trailingIdx {
			lowestIdx = trailingIdx
			lowest = inLow[lowestIdx]
			i = lowestIdx
			i++
			for i <= today {
				tmp = inLow[i]
				if tmp < lowest {
					lowestIdx = i
					lowest = tmp
				}
				i++
			}
			diff = (highest - lowest) / (-100.0)
		} else if tmp <= lowest {
			lowestIdx = today
			lowest = tmp
			diff = (highest - lowest) / (-100.0)
		}
		tmp = inHigh[today]
		if highestIdx < trailingIdx {
			highestIdx = trailingIdx
			highest = inHigh[highestIdx]
			i = highestIdx
			i++
			for i <= today {
				tmp = inHigh[i]
				if tmp > highest {
					highestIdx = i
					highest = tmp
				}
				i++
			}
			diff = (highest - lowest) / (-100.0)
		} else if tmp >= highest {
			highestIdx = today
			highest = tmp
			diff = (highest - lowest) / (-100.0)
		}
		if diff != 0.0 {
			outReal[outIdx] = (highest - inClose[today]) / diff
		} else {
			outReal[outIdx] = 0.0
		}
		outIdx++
		trailingIdx++
		today++
	}
	return outReal
}

/* Volume Indicators */

// Ad - Chaikin A/D Line
func Ad(inHigh []float64, inLow []float64, inClose []float64, inVolume []float64) []float64 {

	outReal := make([]float64, len(inClose))

	startIdx := 0
	nbBar := len(inClose) - startIdx
	currentBar := startIdx
	outIdx := 0
	ad := 0.0
	for nbBar != 0 {
		high := inHigh[currentBar]
		low := inLow[currentBar]
		tmp := high - low
		close := inClose[currentBar]
		if tmp > 0.0 {
			ad += (((close - low) - (high - close)) / tmp) * (inVolume[currentBar])
		}
		outReal[outIdx] = ad
		outIdx++
		currentBar++
		nbBar--
	}
	return outReal
}

// AdOsc - Chaikin A/D Oscillator
func AdOsc(inHigh []float64, inLow []float64, inClose []float64, inVolume []float64, optInFastPeriod int, optInSlowPeriod int) []float64 {

	outReal := make([]float64, len(inClose))

	if (optInFastPeriod < 2) || (optInSlowPeriod < 2) {
		return outReal
	}

	slowestPeriod := 0
	if optInFastPeriod < optInSlowPeriod {
		slowestPeriod = optInSlowPeriod
	} else {
		slowestPeriod = optInFastPeriod
	}
	lookbackTotal := slowestPeriod - 1
	startIdx := lookbackTotal
	today := startIdx - lookbackTotal
	ad := 0.0
	fastk := (2.0 / (float64(optInFastPeriod) + 1.0))
	oneMinusfastk := 1.0 - fastk
	slowk := (2.0 / (float64(optInSlowPeriod) + 1.0))
	oneMinusslowk := 1.0 - slowk
	high := inHigh[today]
	low := inLow[today]
	tmp := high - low
	close := inClose[today]
	if tmp > 0.0 {
		ad += (((close - low) - (high - close)) / tmp) * (inVolume[today])
	}
	today++
	fastEMA := ad
	slowEMA := ad

	for today < startIdx {
		high = inHigh[today]
		low = inLow[today]
		tmp = high - low
		close = inClose[today]
		if tmp > 0.0 {
			ad += (((close - low) - (high - close)) / tmp) * (inVolume[today])
		}
		today++

		fastEMA = (fastk * ad) + (oneMinusfastk * fastEMA)
		slowEMA = (slowk * ad) + (oneMinusslowk * slowEMA)
	}
	outIdx := lookbackTotal
	for today < len(inClose) {
		high = inHigh[today]
		low = inLow[today]
		tmp = high - low
		close = inClose[today]
		if tmp > 0.0 {
			ad += (((close - low) - (high - close)) / tmp) * (inVolume[today])
		}
		today++
		fastEMA = (fastk * ad) + (oneMinusfastk * fastEMA)
		slowEMA = (slowk * ad) + (oneMinusslowk * slowEMA)
		outReal[outIdx] = fastEMA - slowEMA
		outIdx++
	}

	return outReal
}

// Obv - On Balance Volume
func Obv(inReal []float64, inVolume []float64) []float64 {

	outReal := make([]float64, len(inReal))
	startIdx := 0
	prevOBV := inVolume[startIdx]
	prevReal := inReal[startIdx]
	outIdx := 0
	for i := startIdx; i < len(inReal); i++ {
		tempReal := inReal[i]
		if tempReal > prevReal {
			prevOBV += inVolume[i]
		} else if tempReal < prevReal {
			prevOBV -= inVolume[i]
		}
		outReal[outIdx] = prevOBV
		prevReal = tempReal
		outIdx++
	}
	return outReal
}

/* Volatility Indicators */

// Atr - Average True Range
func Atr(inHigh []float64, inLow []float64, inClose []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inClose))

	optInTimePeriodF := float64(optInTimePeriod)

	if optInTimePeriod < 1 {
		return outReal
	}

	if optInTimePeriod <= 1 {
		return TRange(inHigh, inLow, inClose)
	}

	outIdx := optInTimePeriod
	today := optInTimePeriod + 1

	tr := TRange(inHigh, inLow, inClose)
	prevATRTemp := Sma(tr, optInTimePeriod)
	prevATR := prevATRTemp[optInTimePeriod]
	outReal[optInTimePeriod] = prevATR

	for outIdx = optInTimePeriod + 1; outIdx < len(inClose); outIdx++ {
		prevATR *= optInTimePeriodF - 1.0
		prevATR += tr[today]
		prevATR /= optInTimePeriodF
		outReal[outIdx] = prevATR
		today++
	}

	return outReal
}

// Natr - Normalized Average True Range
func Natr(inHigh []float64, inLow []float64, inClose []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inClose))

	if optInTimePeriod < 1 {
		return outReal
	}

	if optInTimePeriod <= 1 {
		return TRange(inHigh, inLow, inClose)
	}

	optInTimePeriodF := float64(optInTimePeriod)
	outIdx := optInTimePeriod
	today := optInTimePeriod

	tr := TRange(inHigh, inLow, inClose)
	prevATRTemp := Sma(tr, optInTimePeriod)
	prevATR := prevATRTemp[optInTimePeriod]

	tempValue := inClose[today]
	if tempValue != 0.0 {
		outReal[outIdx] = (prevATR / tempValue) * 100.0
	} else {
		outReal[outIdx] = 0.0
	}

	for outIdx = optInTimePeriod + 1; outIdx < len(inClose); outIdx++ {
		today++
		prevATR *= optInTimePeriodF - 1.0
		prevATR += tr[today]
		prevATR /= optInTimePeriodF
		tempValue = inClose[today]
		if tempValue != 0.0 {
			outReal[outIdx] = (prevATR / tempValue) * 100.0
		} else {
			outReal[0] = 0.0
		}
	}

	return outReal
}

// TRange - True Range
func TRange(inHigh []float64, inLow []float64, inClose []float64) []float64 {

	outReal := make([]float64, len(inClose))

	startIdx := 1
	outIdx := startIdx
	today := startIdx
	for today < len(inClose) {
		tempLT := inLow[today]
		tempHT := inHigh[today]
		tempCY := inClose[today-1]
		greatest := tempHT - tempLT
		val2 := math.Abs(tempCY - tempHT)
		if val2 > greatest {
			greatest = val2
		}
		val3 := math.Abs(tempCY - tempLT)
		if val3 > greatest {
			greatest = val3
		}
		outReal[outIdx] = greatest
		outIdx++
		today++
	}

	return outReal
}

/* Price Transform */

// AvgPrice - Average Price (o+h+l+c)/4
func AvgPrice(inOpen []float64, inHigh []float64, inLow []float64, inClose []float64) []float64 {

	outReal := make([]float64, len(inClose))
	outIdx := 0
	startIdx := 0

	for i := startIdx; i < len(inClose); i++ {
		outReal[outIdx] = (inHigh[i] + inLow[i] + inClose[i] + inOpen[i]) / 4
		outIdx++
	}
	return outReal
}

// MedPrice - Median Price (h+l)/2
func MedPrice(inHigh []float64, inLow []float64) []float64 {

	outReal := make([]float64, len(inHigh))
	outIdx := 0
	startIdx := 0

	for i := startIdx; i < len(inHigh); i++ {
		outReal[outIdx] = (inHigh[i] + inLow[i]) / 2.0
		outIdx++
	}
	return outReal
}

// TypPrice - Typical Price (h+l+c)/3
func TypPrice(inHigh []float64, inLow []float64, inClose []float64) []float64 {

	outReal := make([]float64, len(inClose))
	outIdx := 0
	startIdx := 0

	for i := startIdx; i < len(inClose); i++ {
		outReal[outIdx] = (inHigh[i] + inLow[i] + inClose[i]) / 3.0
		outIdx++
	}
	return outReal
}

// WclPrice - Weighted Close Price
func WclPrice(inHigh []float64, inLow []float64, inClose []float64) []float64 {

	outReal := make([]float64, len(inClose))
	outIdx := 0
	startIdx := 0

	for i := startIdx; i < len(inClose); i++ {
		outReal[outIdx] = (inHigh[i] + inLow[i] + (inClose[i] * 2.0)) / 4.0
		outIdx++
	}
	return outReal
}

/* Cycle Indicators */

// HtDcPeriod - Hilbert Transform - Dominant Cycle Period (lookback=32)
func HtDcPeriod(inReal []float64) []float64 {

	outReal := make([]float64, len(inReal))

	a := 0.0962
	b := 0.5769
	detrenderOdd := make([]float64, 3)
	detrenderEven := make([]float64, 3)
	q1Odd := make([]float64, 3)
	q1Even := make([]float64, 3)
	jIOdd := make([]float64, 3)
	jIEven := make([]float64, 3)
	jQOdd := make([]float64, 3)
	jQEven := make([]float64, 3)
	rad2Deg := 180.0 / (4.0 * math.Atan(1))
	lookbackTotal := 32
	startIdx := lookbackTotal
	trailingWMAIdx := startIdx - lookbackTotal
	today := trailingWMAIdx
	tempReal := inReal[today]
	today++
	periodWMASub := tempReal
	periodWMASum := tempReal
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 2.0
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 3.0
	trailingWMAValue := 0.0
	i := 9
	smoothedValue := 0.0
	for ok := true; ok; {
		tempReal = inReal[today]
		today++
		periodWMASub += tempReal
		periodWMASub -= trailingWMAValue
		periodWMASum += tempReal * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		smoothedValue = periodWMASum * 0.1
		periodWMASum -= periodWMASub
		i--
		ok = i != 0
	}

	hilbertIdx := 0
	detrender := 0.0
	prevDetrenderOdd := 0.0
	prevDetrenderEven := 0.0
	prevDetrenderInputOdd := 0.0
	prevDetrenderInputEven := 0.0
	Q1 := 0.0
	prevQ1Odd := 0.0
	prevQ1Even := 0.0
	prevQ1InputOdd := 0.0
	prevQ1InputEven := 0.0
	jI := 0.0
	prevJIOdd := 0.0
	prevJIEven := 0.0
	prevJIInputOdd := 0.0
	prevJIInputEven := 0.0
	jQ := 0.0
	prevJQOdd := 0.0
	prevJQEven := 0.0
	prevJQInputOdd := 0.0
	prevJQInputEven := 0.0
	period := 0.0
	outIdx := 32
	prevI2 := 0.0
	prevQ2 := 0.0
	Re := 0.0
	Im := 0.0
	I2 := 0.0
	Q2 := 0.0
	I1ForOddPrev3 := 0.0
	I1ForEvenPrev3 := 0.0
	I1ForOddPrev2 := 0.0
	I1ForEvenPrev2 := 0.0
	smoothPeriod := 0.0
	for today < len(inReal) {
		adjustedPrevPeriod := (0.075 * period) + 0.54
		todayValue := inReal[today]
		periodWMASub += todayValue
		periodWMASub -= trailingWMAValue
		periodWMASum += todayValue * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		smoothedValue = periodWMASum * 0.1
		periodWMASum -= periodWMASub
		hilbertTempReal := 0.0
		if (today % 2) == 0 {
			hilbertTempReal = a * smoothedValue
			detrender = -detrenderEven[hilbertIdx]
			detrenderEven[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderEven
			prevDetrenderEven = b * prevDetrenderInputEven
			detrender += prevDetrenderEven
			prevDetrenderInputEven = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Even[hilbertIdx]
			q1Even[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Even
			prevQ1Even = b * prevQ1InputEven
			Q1 += prevQ1Even
			prevQ1InputEven = detrender
			Q1 *= adjustedPrevPeriod
			hilbertTempReal = a * I1ForEvenPrev3
			jI = -jIEven[hilbertIdx]
			jIEven[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIEven
			prevJIEven = b * prevJIInputEven
			jI += prevJIEven
			prevJIInputEven = I1ForEvenPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQEven[hilbertIdx]
			jQEven[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQEven
			prevJQEven = b * prevJQInputEven
			jQ += prevJQEven
			prevJQInputEven = Q1
			jQ *= adjustedPrevPeriod
			hilbertIdx++
			if hilbertIdx == 3 {
				hilbertIdx = 0
			}
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForEvenPrev3 - jQ)) + (0.8 * prevI2)
			I1ForOddPrev3 = I1ForOddPrev2
			I1ForOddPrev2 = detrender
		} else {
			hilbertTempReal = a * smoothedValue
			detrender = -detrenderOdd[hilbertIdx]
			detrenderOdd[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderOdd
			prevDetrenderOdd = b * prevDetrenderInputOdd
			detrender += prevDetrenderOdd
			prevDetrenderInputOdd = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Odd[hilbertIdx]
			q1Odd[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Odd
			prevQ1Odd = b * prevQ1InputOdd
			Q1 += prevQ1Odd
			prevQ1InputOdd = detrender
			Q1 *= adjustedPrevPeriod
			hilbertTempReal = a * I1ForOddPrev3
			jI = -jIOdd[hilbertIdx]
			jIOdd[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIOdd
			prevJIOdd = b * prevJIInputOdd
			jI += prevJIOdd
			prevJIInputOdd = I1ForOddPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQOdd[hilbertIdx]
			jQOdd[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQOdd
			prevJQOdd = b * prevJQInputOdd
			jQ += prevJQOdd
			prevJQInputOdd = Q1
			jQ *= adjustedPrevPeriod
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForOddPrev3 - jQ)) + (0.8 * prevI2)
			I1ForEvenPrev3 = I1ForEvenPrev2
			I1ForEvenPrev2 = detrender
		}
		Re = (0.2 * ((I2 * prevI2) + (Q2 * prevQ2))) + (0.8 * Re)
		Im = (0.2 * ((I2 * prevQ2) - (Q2 * prevI2))) + (0.8 * Im)
		prevQ2 = Q2
		prevI2 = I2
		tempReal = period
		if (Im != 0.0) && (Re != 0.0) {
			period = 360.0 / (math.Atan(Im/Re) * rad2Deg)
		}
		tempReal2 := 1.5 * tempReal
		if period > tempReal2 {
			period = tempReal2
		}
		tempReal2 = 0.67 * tempReal
		if period < tempReal2 {
			period = tempReal2
		}
		if period < 6 {
			period = 6
		} else if period > 50 {
			period = 50
		}
		period = (0.2 * period) + (0.8 * tempReal)
		smoothPeriod = (0.33 * period) + (0.67 * smoothPeriod)
		if today >= startIdx {
			outReal[outIdx] = smoothPeriod
			outIdx++
		}
		today++
	}
	return outReal
}

// HtDcPhase - Hilbert Transform - Dominant Cycle Phase (lookback=63)
func HtDcPhase(inReal []float64) []float64 {

	outReal := make([]float64, len(inReal))
	a := 0.0962
	b := 0.5769
	detrenderOdd := make([]float64, 3)
	detrenderEven := make([]float64, 3)
	q1Odd := make([]float64, 3)
	q1Even := make([]float64, 3)
	jIOdd := make([]float64, 3)
	jIEven := make([]float64, 3)
	jQOdd := make([]float64, 3)
	jQEven := make([]float64, 3)
	smoothPriceIdx := 0
	maxIdxSmoothPrice := (50 - 1)
	smoothPrice := make([]float64, maxIdxSmoothPrice+1)
	tempReal := math.Atan(1)
	rad2Deg := 45.0 / tempReal
	constDeg2RadBy360 := tempReal * 8.0
	lookbackTotal := 63
	startIdx := lookbackTotal
	trailingWMAIdx := startIdx - lookbackTotal
	today := trailingWMAIdx
	tempReal = inReal[today]
	today++
	periodWMASub := tempReal
	periodWMASum := tempReal
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 2.0
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 3.0
	trailingWMAValue := 0.0
	i := 34
	smoothedValue := 0.0
	for ok := true; ok; {
		tempReal = inReal[today]
		today++
		periodWMASub += tempReal
		periodWMASub -= trailingWMAValue
		periodWMASum += tempReal * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		smoothedValue = periodWMASum * 0.1
		periodWMASum -= periodWMASub
		i--
		ok = i != 0
	}

	hilbertIdx := 0
	detrender := 0.0
	prevDetrenderOdd := 0.0
	prevDetrenderEven := 0.0
	prevDetrenderInputOdd := 0.0
	prevDetrenderInputEven := 0.0
	Q1 := 0.0
	prevQ1Odd := 0.0
	prevQ1Even := 0.0
	prevQ1InputOdd := 0.0
	prevQ1InputEven := 0.0
	jI := 0.0
	prevJIOdd := 0.0
	prevJIEven := 0.0
	prevJIInputOdd := 0.0
	prevJIInputEven := 0.0
	jQ := 0.0
	prevJQOdd := 0.0
	prevJQEven := 0.0
	prevJQInputOdd := 0.0
	prevJQInputEven := 0.0
	period := 0.0
	outIdx := 0
	prevI2 := 0.0
	prevQ2 := 0.0
	Re := 0.0
	Im := 0.0
	I1ForOddPrev3 := 0.0
	I1ForEvenPrev3 := 0.0
	I1ForOddPrev2 := 0.0
	I1ForEvenPrev2 := 0.0
	smoothPeriod := 0.0
	DCPhase := 0.0
	Q2 := 0.0
	I2 := 0.0
	for today < len(inReal) {
		adjustedPrevPeriod := (0.075 * period) + 0.54
		todayValue := inReal[today]
		periodWMASub += todayValue
		periodWMASub -= trailingWMAValue
		periodWMASum += todayValue * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		smoothedValue = periodWMASum * 0.1
		periodWMASum -= periodWMASub
		hilbertTempReal := 0.0
		smoothPrice[smoothPriceIdx] = smoothedValue
		if (today % 2) == 0 {
			hilbertTempReal = a * smoothedValue
			detrender = -detrenderEven[hilbertIdx]
			detrenderEven[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderEven
			prevDetrenderEven = b * prevDetrenderInputEven
			detrender += prevDetrenderEven
			prevDetrenderInputEven = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Even[hilbertIdx]
			q1Even[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Even
			prevQ1Even = b * prevQ1InputEven
			Q1 += prevQ1Even
			prevQ1InputEven = detrender
			Q1 *= adjustedPrevPeriod
			hilbertTempReal = a * I1ForEvenPrev3
			jI = -jIEven[hilbertIdx]
			jIEven[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIEven
			prevJIEven = b * prevJIInputEven
			jI += prevJIEven
			prevJIInputEven = I1ForEvenPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQEven[hilbertIdx]
			jQEven[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQEven
			prevJQEven = b * prevJQInputEven
			jQ += prevJQEven
			prevJQInputEven = Q1
			jQ *= adjustedPrevPeriod
			hilbertIdx++
			if hilbertIdx == 3 {
				hilbertIdx = 0
			}
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForEvenPrev3 - jQ)) + (0.8 * prevI2)
			I1ForOddPrev3 = I1ForOddPrev2
			I1ForOddPrev2 = detrender
		} else {

			hilbertTempReal = a * smoothedValue
			detrender = -detrenderOdd[hilbertIdx]
			detrenderOdd[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderOdd
			prevDetrenderOdd = b * prevDetrenderInputOdd
			detrender += prevDetrenderOdd
			prevDetrenderInputOdd = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Odd[hilbertIdx]
			q1Odd[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Odd
			prevQ1Odd = b * prevQ1InputOdd
			Q1 += prevQ1Odd
			prevQ1InputOdd = detrender
			Q1 *= adjustedPrevPeriod
			hilbertTempReal = a * I1ForOddPrev3
			jI = -jIOdd[hilbertIdx]
			jIOdd[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIOdd
			prevJIOdd = b * prevJIInputOdd
			jI += prevJIOdd
			prevJIInputOdd = I1ForOddPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQOdd[hilbertIdx]
			jQOdd[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQOdd
			prevJQOdd = b * prevJQInputOdd
			jQ += prevJQOdd
			prevJQInputOdd = Q1
			jQ *= adjustedPrevPeriod
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForOddPrev3 - jQ)) + (0.8 * prevI2)
			I1ForEvenPrev3 = I1ForEvenPrev2
			I1ForEvenPrev2 = detrender
		}
		Re = (0.2 * ((I2 * prevI2) + (Q2 * prevQ2))) + (0.8 * Re)
		Im = (0.2 * ((I2 * prevQ2) - (Q2 * prevI2))) + (0.8 * Im)
		prevQ2 = Q2
		prevI2 = I2
		tempReal = period
		if (Im != 0.0) && (Re != 0.0) {
			period = 360.0 / (math.Atan(Im/Re) * rad2Deg)
		}
		tempReal2 := 1.5 * tempReal
		if period > tempReal2 {
			period = tempReal2
		}
		tempReal2 = 0.67 * tempReal
		if period < tempReal2 {
			period = tempReal2
		}
		if period < 6 {
			period = 6
		} else if period > 50 {
			period = 50
		}
		period = (0.2 * period) + (0.8 * tempReal)
		smoothPeriod = (0.33 * period) + (0.67 * smoothPeriod)
		DCPeriod := smoothPeriod + 0.5
		DCPeriodInt := math.Floor(DCPeriod)
		realPart := 0.0
		imagPart := 0.0
		idx := smoothPriceIdx
		for i := 0; i < int(DCPeriodInt); i++ {
			tempReal = (float64(i) * constDeg2RadBy360) / (DCPeriodInt * 1.0)
			tempReal2 = smoothPrice[idx]
			realPart += math.Sin(tempReal) * tempReal2
			imagPart += math.Cos(tempReal) * tempReal2
			if idx == 0 {
				idx = 50 - 1
			} else {
				idx--
			}
		}
		tempReal = math.Abs(imagPart)
		if tempReal > 0.0 {
			DCPhase = math.Atan(realPart/imagPart) * rad2Deg
		} else if tempReal <= 0.01 {
			if realPart < 0.0 {
				DCPhase -= 90.0
			} else if realPart > 0.0 {
				DCPhase += 90.0
			}
		}
		DCPhase += 90.0
		DCPhase += 360.0 / smoothPeriod
		if imagPart < 0.0 {
			DCPhase += 180.0
		}
		if DCPhase > 315.0 {
			DCPhase -= 360.0
		}
		if today >= startIdx {
			outReal[outIdx] = DCPhase
			outIdx++
		}
		smoothPriceIdx++
		if smoothPriceIdx > maxIdxSmoothPrice {
			smoothPriceIdx = 0
		}

		today++
	}
	return outReal
}

// HtPhasor - Hibert Transform - Phasor Components (lookback=32)
func HtPhasor(inReal []float64) ([]float64, []float64) {

	outInPhase := make([]float64, len(inReal))
	outQuadrature := make([]float64, len(inReal))

	a := 0.0962
	b := 0.5769
	detrenderOdd := make([]float64, 3)
	detrenderEven := make([]float64, 3)
	q1Odd := make([]float64, 3)
	q1Even := make([]float64, 3)
	jIOdd := make([]float64, 3)
	jIEven := make([]float64, 3)
	jQOdd := make([]float64, 3)
	jQEven := make([]float64, 3)
	rad2Deg := 180.0 / (4.0 * math.Atan(1))
	lookbackTotal := 32
	startIdx := lookbackTotal
	trailingWMAIdx := startIdx - lookbackTotal
	today := trailingWMAIdx
	tempReal := inReal[today]
	today++
	periodWMASub := tempReal
	periodWMASum := tempReal
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 2.0
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 3.0
	trailingWMAValue := 0.0
	i := 9
	smoothedValue := 0.0
	for ok := true; ok; {
		tempReal = inReal[today]
		today++
		periodWMASub += tempReal
		periodWMASub -= trailingWMAValue
		periodWMASum += tempReal * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		smoothedValue = periodWMASum * 0.1
		periodWMASum -= periodWMASub
		i--
		ok = i != 0
	}
	hilbertIdx := 0
	detrender := 0.0
	prevDetrenderOdd := 0.0
	prevDetrenderEven := 0.0
	prevDetrenderInputOdd := 0.0
	prevDetrenderInputEven := 0.0
	Q1 := 0.0
	prevQ1Odd := 0.0
	prevQ1Even := 0.0
	prevQ1InputOdd := 0.0
	prevQ1InputEven := 0.0
	jI := 0.0
	prevJIOdd := 0.0
	prevJIEven := 0.0
	prevJIInputOdd := 0.0
	prevJIInputEven := 0.0
	jQ := 0.0
	prevJQOdd := 0.0
	prevJQEven := 0.0
	prevJQInputOdd := 0.0
	prevJQInputEven := 0.0
	period := 0.0
	outIdx := 32
	prevI2 := 0.0
	prevQ2 := 0.0
	Re := 0.0
	Im := 0.0
	I1ForOddPrev3 := 0.0
	I1ForEvenPrev3 := 0.0
	I1ForOddPrev2 := 0.0
	I1ForEvenPrev2 := 0.0
	I2 := 0.0
	Q2 := 0.0
	for today < len(inReal) {
		adjustedPrevPeriod := (0.075 * period) + 0.54
		todayValue := inReal[today]
		periodWMASub += todayValue
		periodWMASub -= trailingWMAValue
		periodWMASum += todayValue * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		smoothedValue = periodWMASum * 0.1
		periodWMASum -= periodWMASub
		hilbertTempReal := 0.0
		if (today % 2) == 0 {
			hilbertTempReal = a * smoothedValue
			detrender = -detrenderEven[hilbertIdx]
			detrenderEven[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderEven
			prevDetrenderEven = b * prevDetrenderInputEven
			detrender += prevDetrenderEven
			prevDetrenderInputEven = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Even[hilbertIdx]
			q1Even[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Even
			prevQ1Even = b * prevQ1InputEven
			Q1 += prevQ1Even
			prevQ1InputEven = detrender
			Q1 *= adjustedPrevPeriod

			if today >= startIdx {
				outQuadrature[outIdx] = Q1
				outInPhase[outIdx] = I1ForEvenPrev3
				outIdx++
			}
			hilbertTempReal = a * I1ForEvenPrev3
			jI = -jIEven[hilbertIdx]
			jIEven[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIEven
			prevJIEven = b * prevJIInputEven
			jI += prevJIEven
			prevJIInputEven = I1ForEvenPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQEven[hilbertIdx]
			jQEven[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQEven
			prevJQEven = b * prevJQInputEven
			jQ += prevJQEven
			prevJQInputEven = Q1
			jQ *= adjustedPrevPeriod
			hilbertIdx++
			if hilbertIdx == 3 {
				hilbertIdx = 0
			}
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForEvenPrev3 - jQ)) + (0.8 * prevI2)
			I1ForOddPrev3 = I1ForOddPrev2
			I1ForOddPrev2 = detrender
		} else {

			hilbertTempReal = a * smoothedValue
			detrender = -detrenderOdd[hilbertIdx]
			detrenderOdd[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderOdd
			prevDetrenderOdd = b * prevDetrenderInputOdd
			detrender += prevDetrenderOdd
			prevDetrenderInputOdd = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Odd[hilbertIdx]
			q1Odd[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Odd
			prevQ1Odd = b * prevQ1InputOdd
			Q1 += prevQ1Odd
			prevQ1InputOdd = detrender
			Q1 *= adjustedPrevPeriod
			if today >= startIdx {
				outQuadrature[outIdx] = Q1
				outInPhase[outIdx] = I1ForOddPrev3
				outIdx++
			}
			hilbertTempReal = a * I1ForOddPrev3
			jI = -jIOdd[hilbertIdx]
			jIOdd[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIOdd
			prevJIOdd = b * prevJIInputOdd
			jI += prevJIOdd
			prevJIInputOdd = I1ForOddPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQOdd[hilbertIdx]
			jQOdd[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQOdd
			prevJQOdd = b * prevJQInputOdd
			jQ += prevJQOdd
			prevJQInputOdd = Q1
			jQ *= adjustedPrevPeriod
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForOddPrev3 - jQ)) + (0.8 * prevI2)
			I1ForEvenPrev3 = I1ForEvenPrev2
			I1ForEvenPrev2 = detrender
		}
		Re = (0.2 * ((I2 * prevI2) + (Q2 * prevQ2))) + (0.8 * Re)
		Im = (0.2 * ((I2 * prevQ2) - (Q2 * prevI2))) + (0.8 * Im)
		prevQ2 = Q2
		prevI2 = I2
		tempReal = period
		if (Im != 0.0) && (Re != 0.0) {
			period = 360.0 / (math.Atan(Im/Re) * rad2Deg)
		}
		tempReal2 := 1.5 * tempReal
		if period > tempReal2 {
			period = tempReal2
		}
		tempReal2 = 0.67 * tempReal
		if period < tempReal2 {
			period = tempReal2
		}
		if period < 6 {
			period = 6
		} else if period > 50 {
			period = 50
		}
		period = (0.2 * period) + (0.8 * tempReal)
		today++
	}
	return outInPhase, outQuadrature
}

// HtSine - Hilbert Transform - SineWave (lookback=63)
func HtSine(inReal []float64) ([]float64, []float64) {

	outSine := make([]float64, len(inReal))
	outLeadSine := make([]float64, len(inReal))

	a := 0.0962
	b := 0.5769
	detrenderOdd := make([]float64, 3)
	detrenderEven := make([]float64, 3)
	q1Odd := make([]float64, 3)
	q1Even := make([]float64, 3)
	jIOdd := make([]float64, 3)
	jIEven := make([]float64, 3)
	jQOdd := make([]float64, 3)
	jQEven := make([]float64, 3)
	smoothPriceIdx := 0
	maxIdxSmoothPrice := (50 - 1)
	smoothPrice := make([]float64, maxIdxSmoothPrice+1)
	tempReal := math.Atan(1)
	rad2Deg := 45.0 / tempReal
	deg2Rad := 1.0 / rad2Deg
	constDeg2RadBy360 := tempReal * 8.0
	lookbackTotal := 63
	startIdx := lookbackTotal
	trailingWMAIdx := startIdx - lookbackTotal
	today := trailingWMAIdx
	tempReal = inReal[today]
	today++
	periodWMASub := tempReal
	periodWMASum := tempReal
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 2.0
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 3.0
	trailingWMAValue := 0.0
	i := 34
	smoothedValue := 0.0
	for ok := true; ok; {
		tempReal = inReal[today]
		today++
		periodWMASub += tempReal
		periodWMASub -= trailingWMAValue
		periodWMASum += tempReal * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		smoothedValue = periodWMASum * 0.1
		periodWMASum -= periodWMASub
		i--
		ok = i != 0
	}

	hilbertIdx := 0
	detrender := 0.0
	prevDetrenderOdd := 0.0
	prevDetrenderEven := 0.0
	prevDetrenderInputOdd := 0.0
	prevDetrenderInputEven := 0.0
	Q1 := 0.0
	prevQ1Odd := 0.0
	prevQ1Even := 0.0
	prevQ1InputOdd := 0.0
	prevQ1InputEven := 0.0
	jI := 0.0
	prevJIOdd := 0.0
	prevJIEven := 0.0
	prevJIInputOdd := 0.0
	prevJIInputEven := 0.0
	jQ := 0.0
	prevJQOdd := 0.0
	prevJQEven := 0.0
	prevJQInputOdd := 0.0
	prevJQInputEven := 0.0
	period := 0.0
	outIdx := 63
	prevI2 := 0.0
	prevQ2 := 0.0
	Re := 0.0
	Im := 0.0
	I1ForOddPrev3 := 0.0
	I1ForEvenPrev3 := 0.0
	I1ForOddPrev2 := 0.0
	I1ForEvenPrev2 := 0.0
	smoothPeriod := 0.0
	DCPhase := 0.0
	hilbertTempReal := 0.0
	Q2 := 0.0
	I2 := 0.0
	for today < len(inReal) {
		adjustedPrevPeriod := (0.075 * period) + 0.54
		todayValue := inReal[today]
		periodWMASub += todayValue
		periodWMASub -= trailingWMAValue
		periodWMASum += todayValue * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		smoothedValue = periodWMASum * 0.1
		periodWMASum -= periodWMASub
		smoothPrice[smoothPriceIdx] = smoothedValue
		if (today % 2) == 0 {
			hilbertTempReal = a * smoothedValue
			detrender = -detrenderEven[hilbertIdx]
			detrenderEven[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderEven
			prevDetrenderEven = b * prevDetrenderInputEven
			detrender += prevDetrenderEven
			prevDetrenderInputEven = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Even[hilbertIdx]
			q1Even[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Even
			prevQ1Even = b * prevQ1InputEven
			Q1 += prevQ1Even
			prevQ1InputEven = detrender
			Q1 *= adjustedPrevPeriod
			hilbertTempReal = a * I1ForEvenPrev3
			jI = -jIEven[hilbertIdx]
			jIEven[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIEven
			prevJIEven = b * prevJIInputEven
			jI += prevJIEven
			prevJIInputEven = I1ForEvenPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQEven[hilbertIdx]
			jQEven[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQEven
			prevJQEven = b * prevJQInputEven
			jQ += prevJQEven
			prevJQInputEven = Q1
			jQ *= adjustedPrevPeriod
			hilbertIdx++
			if hilbertIdx == 3 {
				hilbertIdx = 0
			}
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForEvenPrev3 - jQ)) + (0.8 * prevI2)
			I1ForOddPrev3 = I1ForOddPrev2
			I1ForOddPrev2 = detrender
		} else {
			hilbertTempReal = a * smoothedValue
			detrender = -detrenderOdd[hilbertIdx]
			detrenderOdd[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderOdd
			prevDetrenderOdd = b * prevDetrenderInputOdd
			detrender += prevDetrenderOdd
			prevDetrenderInputOdd = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Odd[hilbertIdx]
			q1Odd[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Odd
			prevQ1Odd = b * prevQ1InputOdd
			Q1 += prevQ1Odd
			prevQ1InputOdd = detrender
			Q1 *= adjustedPrevPeriod
			hilbertTempReal = a * I1ForOddPrev3
			jI = -jIOdd[hilbertIdx]
			jIOdd[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIOdd
			prevJIOdd = b * prevJIInputOdd
			jI += prevJIOdd
			prevJIInputOdd = I1ForOddPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQOdd[hilbertIdx]
			jQOdd[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQOdd
			prevJQOdd = b * prevJQInputOdd
			jQ += prevJQOdd
			prevJQInputOdd = Q1
			jQ *= adjustedPrevPeriod
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForOddPrev3 - jQ)) + (0.8 * prevI2)
			I1ForEvenPrev3 = I1ForEvenPrev2
			I1ForEvenPrev2 = detrender
		}
		Re = (0.2 * ((I2 * prevI2) + (Q2 * prevQ2))) + (0.8 * Re)
		Im = (0.2 * ((I2 * prevQ2) - (Q2 * prevI2))) + (0.8 * Im)
		prevQ2 = Q2
		prevI2 = I2
		tempReal = period
		if (Im != 0.0) && (Re != 0.0) {
			period = 360.0 / (math.Atan(Im/Re) * rad2Deg)
		}
		tempReal2 := 1.5 * tempReal
		if period > tempReal2 {
			period = tempReal2
		}
		tempReal2 = 0.67 * tempReal
		if period < tempReal2 {
			period = tempReal2
		}
		if period < 6 {
			period = 6
		} else if period > 50 {
			period = 50
		}
		period = (0.2 * period) + (0.8 * tempReal)
		smoothPeriod = (0.33 * period) + (0.67 * smoothPeriod)
		DCPeriod := smoothPeriod + 0.5
		DCPeriodInt := math.Floor(DCPeriod)
		realPart := 0.0
		imagPart := 0.0
		idx := smoothPriceIdx
		for i := 0; i < int(DCPeriodInt); i++ {
			tempReal = (float64(i) * constDeg2RadBy360) / (DCPeriodInt * 1.0)
			tempReal2 = smoothPrice[idx]
			realPart += math.Sin(tempReal) * tempReal2
			imagPart += math.Cos(tempReal) * tempReal2
			if idx == 0 {
				idx = 50 - 1
			} else {
				idx--
			}
		}
		tempReal = math.Abs(imagPart)
		if tempReal > 0.0 {
			DCPhase = math.Atan(realPart/imagPart) * rad2Deg
		} else if tempReal <= 0.01 {
			if realPart < 0.0 {
				DCPhase -= 90.0
			} else if realPart > 0.0 {
				DCPhase += 90.0
			}
		}
		DCPhase += 90.0
		DCPhase += 360.0 / smoothPeriod
		if imagPart < 0.0 {
			DCPhase += 180.0
		}
		if DCPhase > 315.0 {
			DCPhase -= 360.0
		}
		if today >= startIdx {
			outSine[outIdx] = math.Sin(DCPhase * deg2Rad)
			outLeadSine[outIdx] = math.Sin((DCPhase + 45) * deg2Rad)
			outIdx++
		}
		smoothPriceIdx++
		if smoothPriceIdx > maxIdxSmoothPrice {
			smoothPriceIdx = 0
		}

		today++
	}
	return outSine, outLeadSine
}

// HtTrendline - Hilbert Transform - Trendline (lookback=63)
func HtTrendline(inReal []float64) []float64 {

	outReal := make([]float64, len(inReal))
	a := 0.0962
	b := 0.5769
	detrenderOdd := make([]float64, 3)
	detrenderEven := make([]float64, 3)
	q1Odd := make([]float64, 3)
	q1Even := make([]float64, 3)
	jIOdd := make([]float64, 3)
	jIEven := make([]float64, 3)
	jQOdd := make([]float64, 3)
	jQEven := make([]float64, 3)
	smoothPriceIdx := 0
	maxIdxSmoothPrice := (50 - 1)
	smoothPrice := make([]float64, maxIdxSmoothPrice+1)
	iTrend1 := 0.0
	iTrend2 := 0.0
	iTrend3 := 0.0
	tempReal := math.Atan(1)
	rad2Deg := 45.0 / tempReal
	lookbackTotal := 63
	startIdx := lookbackTotal
	trailingWMAIdx := startIdx - lookbackTotal
	today := trailingWMAIdx
	tempReal = inReal[today]
	today++
	periodWMASub := tempReal
	periodWMASum := tempReal
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 2.0
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 3.0
	trailingWMAValue := 0.0
	i := 34
	for ok := true; ok; {
		tempReal = inReal[today]
		today++
		periodWMASub += tempReal
		periodWMASub -= trailingWMAValue
		periodWMASum += tempReal * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		//smoothedValue := periodWMASum * 0.1
		periodWMASum -= periodWMASub
		i--
		ok = i != 0
	}
	hilbertIdx := 0
	detrender := 0.0
	prevDetrenderOdd := 0.0
	prevDetrenderEven := 0.0
	prevDetrenderInputOdd := 0.0
	prevDetrenderInputEven := 0.0
	Q1 := 0.0
	prevQ1Odd := 0.0
	prevQ1Even := 0.0
	prevQ1InputOdd := 0.0
	prevQ1InputEven := 0.0
	jI := 0.0
	prevJIOdd := 0.0
	prevJIEven := 0.0
	prevJIInputOdd := 0.0
	prevJIInputEven := 0.0
	jQ := 0.0
	prevJQOdd := 0.0
	prevJQEven := 0.0
	prevJQInputOdd := 0.0
	prevJQInputEven := 0.0
	period := 0.0
	outIdx := 63
	prevI2 := 0.0
	prevQ2 := 0.0
	Re := 0.0
	Im := 0.0
	I1ForOddPrev3 := 0.0
	I1ForEvenPrev3 := 0.0
	I1ForOddPrev2 := 0.0
	I1ForEvenPrev2 := 0.0
	smoothPeriod := 0.0
	Q2 := 0.0
	I2 := 0.0
	for today < len(inReal) {
		adjustedPrevPeriod := (0.075 * period) + 0.54
		todayValue := inReal[today]
		periodWMASub += todayValue
		periodWMASub -= trailingWMAValue
		periodWMASum += todayValue * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		smoothedValue := periodWMASum * 0.1
		periodWMASum -= periodWMASub
		smoothPrice[smoothPriceIdx] = smoothedValue
		if (today % 2) == 0 {
			hilbertTempReal := a * smoothedValue
			detrender = -detrenderEven[hilbertIdx]
			detrenderEven[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderEven
			prevDetrenderEven = b * prevDetrenderInputEven
			detrender += prevDetrenderEven
			prevDetrenderInputEven = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Even[hilbertIdx]
			q1Even[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Even
			prevQ1Even = b * prevQ1InputEven
			Q1 += prevQ1Even
			prevQ1InputEven = detrender
			Q1 *= adjustedPrevPeriod
			hilbertTempReal = a * I1ForEvenPrev3
			jI = -jIEven[hilbertIdx]
			jIEven[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIEven
			prevJIEven = b * prevJIInputEven
			jI += prevJIEven
			prevJIInputEven = I1ForEvenPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQEven[hilbertIdx]
			jQEven[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQEven
			prevJQEven = b * prevJQInputEven
			jQ += prevJQEven
			prevJQInputEven = Q1
			jQ *= adjustedPrevPeriod
			hilbertIdx++
			if hilbertIdx == 3 {
				hilbertIdx = 0
			}
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForEvenPrev3 - jQ)) + (0.8 * prevI2)
			I1ForOddPrev3 = I1ForOddPrev2
			I1ForOddPrev2 = detrender
		} else {
			hilbertTempReal := a * smoothedValue
			detrender = -detrenderOdd[hilbertIdx]
			detrenderOdd[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderOdd
			prevDetrenderOdd = b * prevDetrenderInputOdd
			detrender += prevDetrenderOdd
			prevDetrenderInputOdd = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Odd[hilbertIdx]
			q1Odd[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Odd
			prevQ1Odd = b * prevQ1InputOdd
			Q1 += prevQ1Odd
			prevQ1InputOdd = detrender
			Q1 *= adjustedPrevPeriod
			hilbertTempReal = a * I1ForOddPrev3
			jI = -jIOdd[hilbertIdx]
			jIOdd[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIOdd
			prevJIOdd = b * prevJIInputOdd
			jI += prevJIOdd
			prevJIInputOdd = I1ForOddPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQOdd[hilbertIdx]
			jQOdd[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQOdd
			prevJQOdd = b * prevJQInputOdd
			jQ += prevJQOdd
			prevJQInputOdd = Q1
			jQ *= adjustedPrevPeriod
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForOddPrev3 - jQ)) + (0.8 * prevI2)
			I1ForEvenPrev3 = I1ForEvenPrev2
			I1ForEvenPrev2 = detrender
		}
		Re = (0.2 * ((I2 * prevI2) + (Q2 * prevQ2))) + (0.8 * Re)
		Im = (0.2 * ((I2 * prevQ2) - (Q2 * prevI2))) + (0.8 * Im)
		prevQ2 = Q2
		prevI2 = I2
		tempReal = period
		if (Im != 0.0) && (Re != 0.0) {
			period = 360.0 / (math.Atan(Im/Re) * rad2Deg)
		}
		tempReal2 := 1.5 * tempReal
		if period > tempReal2 {
			period = tempReal2
		}
		tempReal2 = 0.67 * tempReal
		if period < tempReal2 {
			period = tempReal2
		}
		if period < 6 {
			period = 6
		} else if period > 50 {
			period = 50
		}
		period = (0.2 * period) + (0.8 * tempReal)
		smoothPeriod = (0.33 * period) + (0.67 * smoothPeriod)
		DCPeriod := smoothPeriod + 0.5
		DCPeriodInt := math.Floor(DCPeriod)
		idx := today
		tempReal = 0.0
		for i := 0; i < int(DCPeriodInt); i++ {
			tempReal += inReal[idx]
			idx--
		}
		if DCPeriodInt > 0 {
			tempReal = tempReal / (DCPeriodInt * 1.0)
		}
		tempReal2 = (4.0*tempReal + 3.0*iTrend1 + 2.0*iTrend2 + iTrend3) / 10.0
		iTrend3 = iTrend2
		iTrend2 = iTrend1
		iTrend1 = tempReal
		if today >= startIdx {
			outReal[outIdx] = tempReal2
			outIdx++
		}
		smoothPriceIdx++
		if smoothPriceIdx > maxIdxSmoothPrice {
			smoothPriceIdx = 0
		}

		today++
	}
	return outReal
}

// HtTrendMode - Hilbert Transform - Trend vs Cycle Mode (lookback=63)
func HtTrendMode(inReal []float64) []float64 {

	outReal := make([]float64, len(inReal))
	a := 0.0962
	b := 0.5769
	detrenderOdd := make([]float64, 3)
	detrenderEven := make([]float64, 3)
	q1Odd := make([]float64, 3)
	q1Even := make([]float64, 3)
	jIOdd := make([]float64, 3)
	jIEven := make([]float64, 3)
	jQOdd := make([]float64, 3)
	jQEven := make([]float64, 3)
	smoothPriceIdx := 0
	maxIdxSmoothPrice := (50 - 1)
	smoothPrice := make([]float64, maxIdxSmoothPrice+1)
	iTrend1 := 0.0
	iTrend2 := 0.0
	iTrend3 := 0.0
	daysInTrend := 0
	prevDCPhase := 0.0
	DCPhase := 0.0
	prevSine := 0.0
	sine := 0.0
	prevLeadSine := 0.0
	leadSine := 0.0
	tempReal := math.Atan(1)
	rad2Deg := 45.0 / tempReal
	deg2Rad := 1.0 / rad2Deg
	constDeg2RadBy360 := tempReal * 8.0
	lookbackTotal := 63
	startIdx := lookbackTotal
	trailingWMAIdx := startIdx - lookbackTotal
	today := trailingWMAIdx
	tempReal = inReal[today]
	today++
	periodWMASub := tempReal
	periodWMASum := tempReal
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 2.0
	tempReal = inReal[today]
	today++
	periodWMASub += tempReal
	periodWMASum += tempReal * 3.0
	trailingWMAValue := 0.0
	i := 34

	for ok := true; ok; {
		tempReal = inReal[today]
		today++
		periodWMASub += tempReal
		periodWMASub -= trailingWMAValue
		periodWMASum += tempReal * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		//smoothedValue := periodWMASum * 0.1
		periodWMASum -= periodWMASub
		i--
		ok = i != 0
	}

	hilbertIdx := 0
	detrender := 0.0
	prevDetrenderOdd := 0.0
	prevDetrenderEven := 0.0
	prevDetrenderInputOdd := 0.0
	prevDetrenderInputEven := 0.0
	Q1 := 0.0
	prevQ1Odd := 0.0
	prevQ1Even := 0.0
	prevQ1InputOdd := 0.0
	prevQ1InputEven := 0.0
	jI := 0.0
	prevJIOdd := 0.0
	prevJIEven := 0.0
	prevJIInputOdd := 0.0
	prevJIInputEven := 0.0
	jQ := 0.0
	prevJQOdd := 0.0
	prevJQEven := 0.0
	prevJQInputOdd := 0.0
	prevJQInputEven := 0.0
	period := 0.0
	outIdx := 63
	prevI2 := 0.0
	prevQ2 := 0.0
	Re := 0.0
	Im := 0.0
	I1ForOddPrev3 := 0.0
	I1ForEvenPrev3 := 0.0
	I1ForOddPrev2 := 0.0
	I1ForEvenPrev2 := 0.0
	smoothPeriod := 0.0
	DCPhase = 0.0
	smoothedValue := 0.0
	hilbertTempReal := 0.0
	Q2 := 0.0
	I2 := 0.0
	for today < len(inReal) {
		adjustedPrevPeriod := (0.075 * period) + 0.54
		todayValue := inReal[today]
		periodWMASub += todayValue
		periodWMASub -= trailingWMAValue
		periodWMASum += todayValue * 4.0
		trailingWMAValue = inReal[trailingWMAIdx]
		trailingWMAIdx++
		smoothedValue = periodWMASum * 0.1
		periodWMASum -= periodWMASub

		smoothPrice[smoothPriceIdx] = smoothedValue
		if (today % 2) == 0 {
			hilbertTempReal = a * smoothedValue
			detrender = -detrenderEven[hilbertIdx]
			detrenderEven[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderEven
			prevDetrenderEven = b * prevDetrenderInputEven
			detrender += prevDetrenderEven
			prevDetrenderInputEven = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Even[hilbertIdx]
			q1Even[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Even
			prevQ1Even = b * prevQ1InputEven
			Q1 += prevQ1Even
			prevQ1InputEven = detrender
			Q1 *= adjustedPrevPeriod
			hilbertTempReal = a * I1ForEvenPrev3
			jI = -jIEven[hilbertIdx]
			jIEven[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIEven
			prevJIEven = b * prevJIInputEven
			jI += prevJIEven
			prevJIInputEven = I1ForEvenPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQEven[hilbertIdx]
			jQEven[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQEven
			prevJQEven = b * prevJQInputEven
			jQ += prevJQEven
			prevJQInputEven = Q1
			jQ *= adjustedPrevPeriod
			hilbertIdx++
			if hilbertIdx == 3 {
				hilbertIdx = 0
			}
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForEvenPrev3 - jQ)) + (0.8 * prevI2)
			I1ForOddPrev3 = I1ForOddPrev2
			I1ForOddPrev2 = detrender
		} else {
			hilbertTempReal = a * smoothedValue
			detrender = -detrenderOdd[hilbertIdx]
			detrenderOdd[hilbertIdx] = hilbertTempReal
			detrender += hilbertTempReal
			detrender -= prevDetrenderOdd
			prevDetrenderOdd = b * prevDetrenderInputOdd
			detrender += prevDetrenderOdd
			prevDetrenderInputOdd = smoothedValue
			detrender *= adjustedPrevPeriod
			hilbertTempReal = a * detrender
			Q1 = -q1Odd[hilbertIdx]
			q1Odd[hilbertIdx] = hilbertTempReal
			Q1 += hilbertTempReal
			Q1 -= prevQ1Odd
			prevQ1Odd = b * prevQ1InputOdd
			Q1 += prevQ1Odd
			prevQ1InputOdd = detrender
			Q1 *= adjustedPrevPeriod
			hilbertTempReal = a * I1ForOddPrev3
			jI = -jIOdd[hilbertIdx]
			jIOdd[hilbertIdx] = hilbertTempReal
			jI += hilbertTempReal
			jI -= prevJIOdd
			prevJIOdd = b * prevJIInputOdd
			jI += prevJIOdd
			prevJIInputOdd = I1ForOddPrev3
			jI *= adjustedPrevPeriod
			hilbertTempReal = a * Q1
			jQ = -jQOdd[hilbertIdx]
			jQOdd[hilbertIdx] = hilbertTempReal
			jQ += hilbertTempReal
			jQ -= prevJQOdd
			prevJQOdd = b * prevJQInputOdd
			jQ += prevJQOdd
			prevJQInputOdd = Q1
			jQ *= adjustedPrevPeriod
			Q2 = (0.2 * (Q1 + jI)) + (0.8 * prevQ2)
			I2 = (0.2 * (I1ForOddPrev3 - jQ)) + (0.8 * prevI2)
			I1ForEvenPrev3 = I1ForEvenPrev2
			I1ForEvenPrev2 = detrender
		}
		Re = (0.2 * ((I2 * prevI2) + (Q2 * prevQ2))) + (0.8 * Re)
		Im = (0.2 * ((I2 * prevQ2) - (Q2 * prevI2))) + (0.8 * Im)
		prevQ2 = Q2
		prevI2 = I2
		tempReal = period
		if (Im != 0.0) && (Re != 0.0) {
			period = 360.0 / (math.Atan(Im/Re) * rad2Deg)
		}
		tempReal2 := 1.5 * tempReal
		if period > tempReal2 {
			period = tempReal2
		}
		tempReal2 = 0.67 * tempReal
		if period < tempReal2 {
			period = tempReal2
		}
		if period < 6 {
			period = 6
		} else if period > 50 {
			period = 50
		}
		period = (0.2 * period) + (0.8 * tempReal)
		smoothPeriod = (0.33 * period) + (0.67 * smoothPeriod)
		prevDCPhase = DCPhase
		DCPeriod := smoothPeriod + 0.5
		DCPeriodInt := math.Floor(DCPeriod)
		realPart := 0.0
		imagPart := 0.0
		idx := smoothPriceIdx
		for i := 0; i < int(DCPeriodInt); i++ {
			tempReal = (float64(i) * constDeg2RadBy360) / (DCPeriodInt * 1.0)
			tempReal2 = smoothPrice[idx]
			realPart += math.Sin(tempReal) * tempReal2
			imagPart += math.Cos(tempReal) * tempReal2
			if idx == 0 {
				idx = 50 - 1
			} else {
				idx--
			}
		}
		tempReal = math.Abs(imagPart)
		if tempReal > 0.0 {
			DCPhase = math.Atan(realPart/imagPart) * rad2Deg
		} else if tempReal <= 0.01 {
			if realPart < 0.0 {
				DCPhase -= 90.0
			} else if realPart > 0.0 {
				DCPhase += 90.0
			}
		}
		DCPhase += 90.0
		DCPhase += 360.0 / smoothPeriod
		if imagPart < 0.0 {
			DCPhase += 180.0
		}
		if DCPhase > 315.0 {
			DCPhase -= 360.0
		}
		prevSine = sine
		prevLeadSine = leadSine
		sine = math.Sin(DCPhase * deg2Rad)
		leadSine = math.Sin((DCPhase + 45) * deg2Rad)
		DCPeriod = smoothPeriod + 0.5
		DCPeriodInt = math.Floor(DCPeriod)
		idx = today
		tempReal = 0.0
		for i := 0; i < int(DCPeriodInt); i++ {
			tempReal += inReal[idx]
			idx--
		}
		if DCPeriodInt > 0 {
			tempReal = tempReal / (DCPeriodInt * 1.0)
		}
		trendline := (4.0*tempReal + 3.0*iTrend1 + 2.0*iTrend2 + iTrend3) / 10.0
		iTrend3 = iTrend2
		iTrend2 = iTrend1
		iTrend1 = tempReal
		trend := 1
		if ((sine > leadSine) && (prevSine <= prevLeadSine)) || ((sine < leadSine) && (prevSine >= prevLeadSine)) {
			daysInTrend = 0
			trend = 0
		}
		daysInTrend++
		if float64(daysInTrend) < (0.5 * smoothPeriod) {
			trend = 0
		}
		tempReal = DCPhase - prevDCPhase
		if (smoothPeriod != 0.0) && ((tempReal > (0.67 * 360.0 / smoothPeriod)) && (tempReal < (1.5 * 360.0 / smoothPeriod))) {
			trend = 0
		}
		tempReal = smoothPrice[smoothPriceIdx]
		if (trendline != 0.0) && (math.Abs((tempReal-trendline)/trendline) >= 0.015) {
			trend = 1
		}
		if today >= startIdx {
			outReal[outIdx] = float64(trend)
			outIdx++
		}
		smoothPriceIdx++
		if smoothPriceIdx > maxIdxSmoothPrice {
			smoothPriceIdx = 0
		}

		today++
	}
	return outReal
}

/* Statistic Functions */

// Beta - Beta
func Beta(inReal0 []float64, inReal1 []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal0))

	x := 0.0
	y := 0.0
	sSS := 0.0
	sXY := 0.0
	sX := 0.0
	sY := 0.0
	tmpReal := 0.0
	n := 0.0
	nbInitialElementNeeded := optInTimePeriod
	startIdx := nbInitialElementNeeded
	trailingIdx := startIdx - nbInitialElementNeeded
	trailingLastPriceX := inReal0[trailingIdx]
	lastPriceX := trailingLastPriceX
	trailingLastPriceY := inReal1[trailingIdx]
	lastPriceY := trailingLastPriceY
	trailingIdx++
	i := trailingIdx
	for i < startIdx {
		tmpReal := inReal0[i]
		x := 0.0
		if !((-0.00000000000001 < lastPriceX) && (lastPriceX < 0.00000000000001)) {
			x = (tmpReal - lastPriceX) / lastPriceX
		}
		lastPriceX = tmpReal
		tmpReal = inReal1[i]
		i++
		y := 0.0
		if !((-0.00000000000001 < lastPriceY) && (lastPriceY < 0.00000000000001)) {
			y = (tmpReal - lastPriceY) / lastPriceY
		}
		lastPriceY = tmpReal
		sSS += x * x
		sXY += x * y
		sX += x
		sY += y
	}
	outIdx := optInTimePeriod
	n = float64(optInTimePeriod)
	for ok := true; ok; {
		tmpReal = inReal0[i]
		if !((-0.00000000000001 < lastPriceX) && (lastPriceX < 0.00000000000001)) {
			x = (tmpReal - lastPriceX) / lastPriceX
		} else {
			x = 0.0
		}
		lastPriceX = tmpReal
		tmpReal = inReal1[i]
		i++
		if !((-0.00000000000001 < lastPriceY) && (lastPriceY < 0.00000000000001)) {
			y = (tmpReal - lastPriceY) / lastPriceY
		} else {
			y = 0.0
		}
		lastPriceY = tmpReal
		sSS += x * x
		sXY += x * y
		sX += x
		sY += y
		tmpReal = inReal0[trailingIdx]
		if !(((-(0.00000000000001)) < trailingLastPriceX) && (trailingLastPriceX < (0.00000000000001))) {
			x = (tmpReal - trailingLastPriceX) / trailingLastPriceX
		} else {
			x = 0.0
		}
		trailingLastPriceX = tmpReal
		tmpReal = inReal1[trailingIdx]
		trailingIdx++
		if !(((-(0.00000000000001)) < trailingLastPriceY) && (trailingLastPriceY < (0.00000000000001))) {
			y = (tmpReal - trailingLastPriceY) / trailingLastPriceY
		} else {
			y = 0.0
		}
		trailingLastPriceY = tmpReal
		tmpReal = (n * sSS) - (sX * sX)
		if !(((-(0.00000000000001)) < tmpReal) && (tmpReal < (0.00000000000001))) {
			outReal[outIdx] = ((n * sXY) - (sX * sY)) / tmpReal
		} else {
			outReal[outIdx] = 0.0
		}
		outIdx++
		sSS -= x * x
		sXY -= x * y
		sX -= x
		sY -= y
		ok = i < len(inReal0)
	}

	return outReal
}

// Correl - Pearson's Correlation Coefficient (r)
func Correl(inReal0 []float64, inReal1 []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal0))

	optInTimePeriodF := float64(optInTimePeriod)
	lookbackTotal := optInTimePeriod - 1
	startIdx := lookbackTotal
	trailingIdx := startIdx - lookbackTotal
	sumXY, sumX, sumY, sumX2, sumY2 := 0.0, 0.0, 0.0, 0.0, 0.0
	today := trailingIdx
	for today = trailingIdx; today <= startIdx; today++ {
		x := inReal0[today]
		sumX += x
		sumX2 += x * x
		y := inReal1[today]
		sumXY += x * y
		sumY += y
		sumY2 += y * y
	}
	trailingX := inReal0[trailingIdx]
	trailingY := inReal1[trailingIdx]
	trailingIdx++
	tempReal := (sumX2 - ((sumX * sumX) / optInTimePeriodF)) * (sumY2 - ((sumY * sumY) / optInTimePeriodF))
	if !(tempReal < 0.00000000000001) {
		outReal[optInTimePeriod-1] = (sumXY - ((sumX * sumY) / optInTimePeriodF)) / math.Sqrt(tempReal)
	} else {
		outReal[optInTimePeriod-1] = 0.0
	}
	outIdx := optInTimePeriod
	for today < len(inReal0) {
		sumX -= trailingX
		sumX2 -= trailingX * trailingX
		sumXY -= trailingX * trailingY
		sumY -= trailingY
		sumY2 -= trailingY * trailingY
		x := inReal0[today]
		sumX += x
		sumX2 += x * x
		y := inReal1[today]
		today++
		sumXY += x * y
		sumY += y
		sumY2 += y * y
		trailingX = inReal0[trailingIdx]
		trailingY = inReal1[trailingIdx]
		trailingIdx++
		tempReal = (sumX2 - ((sumX * sumX) / optInTimePeriodF)) * (sumY2 - ((sumY * sumY) / optInTimePeriodF))
		if !(tempReal < (0.00000000000001)) {
			outReal[outIdx] = (sumXY - ((sumX * sumY) / optInTimePeriodF)) / math.Sqrt(tempReal)
		} else {
			outReal[outIdx] = 0.0
		}
		outIdx++
	}
	return outReal
}

// LinearReg - Linear Regression
func LinearReg(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	optInTimePeriodF := float64(optInTimePeriod)
	lookbackTotal := optInTimePeriod
	startIdx := lookbackTotal
	outIdx := startIdx - 1
	today := startIdx - 1
	sumX := optInTimePeriodF * (optInTimePeriodF - 1) * 0.5
	sumXSqr := optInTimePeriodF * (optInTimePeriodF - 1) * (2*optInTimePeriodF - 1) / 6
	divisor := sumX*sumX - optInTimePeriodF*sumXSqr
	for today < len(inReal) {
		sumXY := 0.0
		sumY := 0.0
		i := optInTimePeriod
		for i != 0 {
			i--
			tempValue1 := inReal[today-i]
			sumY += tempValue1
			sumXY += float64(i) * tempValue1
		}
		m := (optInTimePeriodF*sumXY - sumX*sumY) / divisor
		b := (sumY - m*sumX) / optInTimePeriodF
		outReal[outIdx] = b + m*(optInTimePeriodF-1)
		outIdx++
		today++
	}
	return outReal
}

// LinearRegAngle - Linear Regression Angle
func LinearRegAngle(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	optInTimePeriodF := float64(optInTimePeriod)
	lookbackTotal := optInTimePeriod
	startIdx := lookbackTotal
	outIdx := startIdx - 1
	today := startIdx - 1
	sumX := optInTimePeriodF * (optInTimePeriodF - 1) * 0.5
	sumXSqr := optInTimePeriodF * (optInTimePeriodF - 1) * (2*optInTimePeriodF - 1) / 6
	divisor := sumX*sumX - optInTimePeriodF*sumXSqr
	for today < len(inReal) {
		sumXY := 0.0
		sumY := 0.0
		i := optInTimePeriod
		for i != 0 {
			i--
			tempValue1 := inReal[today-i]
			sumY += tempValue1
			sumXY += float64(i) * tempValue1
		}
		m := (optInTimePeriodF*sumXY - sumX*sumY) / divisor
		outReal[outIdx] = math.Atan(m) * (180.0 / 3.14159265358979323846)
		outIdx++
		today++
	}
	return outReal
}

// LinearRegIntercept - Linear Regression Intercept
func LinearRegIntercept(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	optInTimePeriodF := float64(optInTimePeriod)
	lookbackTotal := optInTimePeriod
	startIdx := lookbackTotal
	outIdx := startIdx - 1
	today := startIdx - 1
	sumX := optInTimePeriodF * (optInTimePeriodF - 1) * 0.5
	sumXSqr := optInTimePeriodF * (optInTimePeriodF - 1) * (2*optInTimePeriodF - 1) / 6
	divisor := sumX*sumX - optInTimePeriodF*sumXSqr
	for today < len(inReal) {
		sumXY := 0.0
		sumY := 0.0
		i := optInTimePeriod
		for i != 0 {
			i--
			tempValue1 := inReal[today-i]
			sumY += tempValue1
			sumXY += float64(i) * tempValue1
		}
		m := (optInTimePeriodF*sumXY - sumX*sumY) / divisor
		outReal[outIdx] = (sumY - m*sumX) / optInTimePeriodF
		outIdx++
		today++
	}
	return outReal
}

// LinearRegSlope - Linear Regression Slope
func LinearRegSlope(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	optInTimePeriodF := float64(optInTimePeriod)
	lookbackTotal := optInTimePeriod
	startIdx := lookbackTotal
	outIdx := startIdx - 1
	today := startIdx - 1
	sumX := optInTimePeriodF * (optInTimePeriodF - 1) * 0.5
	sumXSqr := optInTimePeriodF * (optInTimePeriodF - 1) * (2*optInTimePeriodF - 1) / 6
	divisor := sumX*sumX - optInTimePeriodF*sumXSqr
	for today < len(inReal) {
		sumXY := 0.0
		sumY := 0.0
		i := optInTimePeriod
		for i != 0 {
			i--
			tempValue1 := inReal[today-i]
			sumY += tempValue1
			sumXY += float64(i) * tempValue1
		}
		outReal[outIdx] = (optInTimePeriodF*sumXY - sumX*sumY) / divisor
		outIdx++
		today++
	}
	return outReal
}

// StdDev - Standard Deviation
func StdDev(inReal []float64, optInTimePeriod int, optInNbDev float64) []float64 {

	outReal := Var(inReal, optInTimePeriod)

	if optInNbDev != 1.0 {
		for i := 0; i < len(inReal); i++ {
			tempReal := outReal[i]
			if !(tempReal < 0.00000000000001) {
				outReal[i] = math.Sqrt(tempReal) * optInNbDev
			} else {
				outReal[i] = 0.0
			}
		}
	} else {
		for i := 0; i < len(inReal); i++ {
			tempReal := outReal[i]
			if !(tempReal < 0.00000000000001) {
				outReal[i] = math.Sqrt(tempReal)
			} else {
				outReal[i] = 0.0
			}
		}
	}
	return outReal
}

// Tsf - Time Series Forecast
func Tsf(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	optInTimePeriodF := float64(optInTimePeriod)
	lookbackTotal := optInTimePeriod
	startIdx := lookbackTotal
	outIdx := startIdx - 1
	today := startIdx - 1
	sumX := optInTimePeriodF * (optInTimePeriodF - 1.0) * 0.5
	sumXSqr := optInTimePeriodF * (optInTimePeriodF - 1) * (2*optInTimePeriodF - 1) / 6
	divisor := sumX*sumX - optInTimePeriodF*sumXSqr
	for today < len(inReal) {
		sumXY := 0.0
		sumY := 0.0
		i := optInTimePeriod
		for i != 0 {
			i--
			tempValue1 := inReal[today-i]
			sumY += tempValue1
			sumXY += float64(i) * tempValue1
		}
		m := (optInTimePeriodF*sumXY - sumX*sumY) / divisor
		b := (sumY - m*sumX) / optInTimePeriodF
		outReal[outIdx] = b + m*optInTimePeriodF
		today++
		outIdx++
	}
	return outReal
}

// Var - Variance
func Var(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	nbInitialElementNeeded := optInTimePeriod - 1
	startIdx := nbInitialElementNeeded
	periodTotal1 := 0.0
	periodTotal2 := 0.0
	trailingIdx := startIdx - nbInitialElementNeeded
	i := trailingIdx
	if optInTimePeriod > 1 {
		for i < startIdx {
			tempReal := inReal[i]
			periodTotal1 += tempReal
			tempReal *= tempReal
			periodTotal2 += tempReal
			i++
		}
	}
	outIdx := startIdx
	for ok := true; ok; {
		tempReal := inReal[i]
		periodTotal1 += tempReal
		tempReal *= tempReal
		periodTotal2 += tempReal
		meanValue1 := periodTotal1 / float64(optInTimePeriod)
		meanValue2 := periodTotal2 / float64(optInTimePeriod)
		tempReal = inReal[trailingIdx]
		periodTotal1 -= tempReal
		tempReal *= tempReal
		periodTotal2 -= tempReal
		outReal[outIdx] = meanValue2 - meanValue1*meanValue1
		i++
		trailingIdx++
		outIdx++
		ok = i < len(inReal)
	}
	return outReal
}

/* Math Transform Functions */

// Acos - Vector Trigonometric ACOS
func Acos(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Acos(inReal[i])
	}
	return outReal
}

// Asin - Vector Trigonometric ASIN
func Asin(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Asin(inReal[i])
	}
	return outReal
}

// Atan - Vector Trigonometric ATAN
func Atan(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Atan(inReal[i])
	}
	return outReal
}

// Ceil - Vector CEIL
func Ceil(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Ceil(inReal[i])
	}
	return outReal
}

// Cos - Vector Trigonometric COS
func Cos(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Cos(inReal[i])
	}
	return outReal
}

// Cosh - Vector Trigonometric COSH
func Cosh(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Cosh(inReal[i])
	}
	return outReal
}

// Exp - Vector atrithmetic EXP
func Exp(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Exp(inReal[i])
	}
	return outReal
}

// Floor - Vector FLOOR
func Floor(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Floor(inReal[i])
	}
	return outReal
}

// Ln - Vector natural log LN
func Ln(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Log(inReal[i])
	}
	return outReal
}

// Log10 - Vector LOG10
func Log10(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Log10(inReal[i])
	}
	return outReal
}

// Sin - Vector Trigonometric SIN
func Sin(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Sin(inReal[i])
	}
	return outReal
}

// Sinh - Vector Trigonometric SINH
func Sinh(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Sinh(inReal[i])
	}
	return outReal
}

// Sqrt - Vector SQRT
func Sqrt(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Sqrt(inReal[i])
	}
	return outReal
}

// Tan - Vector Trigonometric TAN
func Tan(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Tan(inReal[i])
	}
	return outReal
}

// Tanh - Vector Trigonometric TANH
func Tanh(inReal []float64) []float64 {
	outReal := make([]float64, len(inReal))
	for i := 0; i < len(inReal); i++ {
		outReal[i] = math.Tanh(inReal[i])
	}
	return outReal
}

/* Math Operator Functions
TODO:
  MINMAX - Lowest and highest values over a specified period
    min, max = MINMAX(close, timeperiod=30)
  MINMAXINDEX - Indexes of lowest and highest values over a specified period
    minidx, maxidx = MINMAXINDEX(close, timeperiod=30)
*/

// Add - Vector arithmetic addition
func Add(inReal0 []float64, inReal1 []float64) []float64 {
	outReal := make([]float64, len(inReal0))
	for i := 0; i < len(inReal0); i++ {
		outReal[i] = inReal0[i] + inReal1[i]
	}
	return outReal
}

// Div - Vector arithmetic division
func Div(inReal0 []float64, inReal1 []float64) []float64 {
	outReal := make([]float64, len(inReal0))
	for i := 0; i < len(inReal0); i++ {
		outReal[i] = inReal0[i] / inReal1[i]
	}
	return outReal
}

// Max - Highest value over a period
func Max(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	if optInTimePeriod < 2 {
		return outReal
	}

	nbInitialElementNeeded := optInTimePeriod - 1
	startIdx := nbInitialElementNeeded
	outIdx := startIdx
	today := startIdx
	trailingIdx := startIdx - nbInitialElementNeeded
	highestIdx := -1
	highest := 0.0

	for today < len(outReal) {

		tmp := inReal[today]

		if highestIdx < trailingIdx {
			highestIdx = trailingIdx
			highest = inReal[highestIdx]
			i := highestIdx + 1
			for i <= today {
				tmp = inReal[i]
				if tmp > highest {
					highestIdx = i
					highest = tmp
				}
				i++
			}
		} else if tmp >= highest {
			highestIdx = today
			highest = tmp
		}
		outReal[outIdx] = highest
		outIdx++
		trailingIdx++
		today++
	}

	return outReal
}

// MaxIndex - Index of highest value over a specified period
func MaxIndex(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	if optInTimePeriod < 2 {
		return outReal
	}

	nbInitialElementNeeded := optInTimePeriod - 1
	startIdx := nbInitialElementNeeded
	outIdx := startIdx
	today := startIdx
	trailingIdx := startIdx - nbInitialElementNeeded
	highestIdx := -1
	highest := 0.0
	for today < len(inReal) {
		tmp := inReal[today]
		if highestIdx < trailingIdx {
			highestIdx = trailingIdx
			highest = inReal[highestIdx]
			i := highestIdx + 1
			for i <= today {
				tmp := inReal[i]
				if tmp > highest {
					highestIdx = i
					highest = tmp
				}
				i++
			}
		} else if tmp >= highest {
			highestIdx = today
			highest = tmp
		}
		outReal[outIdx] = float64(highestIdx)
		outIdx++
		trailingIdx++
		today++
	}

	return outReal
}

// Min - Lowest value over a period
func Min(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	if optInTimePeriod < 2 {
		return outReal
	}

	nbInitialElementNeeded := optInTimePeriod - 1
	startIdx := nbInitialElementNeeded
	outIdx := startIdx
	today := startIdx
	trailingIdx := startIdx - nbInitialElementNeeded
	lowestIdx := -1
	lowest := 0.0
	for today < len(outReal) {

		tmp := inReal[today]

		if lowestIdx < trailingIdx {
			lowestIdx = trailingIdx
			lowest = inReal[lowestIdx]
			i := lowestIdx + 1
			for i <= today {
				tmp = inReal[i]
				if tmp < lowest {
					lowestIdx = i
					lowest = tmp
				}
				i++
			}
		} else if tmp <= lowest {
			lowestIdx = today
			lowest = tmp
		}
		outReal[outIdx] = lowest
		outIdx++
		trailingIdx++
		today++
	}

	return outReal
}

// MinIndex - Index of lowest value over a specified period
func MinIndex(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	if optInTimePeriod < 2 {
		return outReal
	}

	nbInitialElementNeeded := optInTimePeriod - 1
	startIdx := nbInitialElementNeeded
	outIdx := startIdx
	today := startIdx
	trailingIdx := startIdx - nbInitialElementNeeded
	lowestIdx := -1
	lowest := 0.0
	for today < len(inReal) {
		tmp := inReal[today]
		if lowestIdx < trailingIdx {
			lowestIdx = trailingIdx
			lowest = inReal[lowestIdx]
			i := lowestIdx + 1
			for i <= today {
				tmp = inReal[i]
				if tmp < lowest {
					lowestIdx = i
					lowest = tmp
				}
				i++
			}
		} else if tmp <= lowest {
			lowestIdx = today
			lowest = tmp
		}
		outReal[outIdx] = float64(lowestIdx)
		outIdx++
		trailingIdx++
		today++
	}
	return outReal
}

// Mult - Vector arithmetic multiply
func Mult(inReal0 []float64, inReal1 []float64) []float64 {
	outReal := make([]float64, len(inReal0))
	for i := 0; i < len(inReal0); i++ {
		outReal[i] = inReal0[i] * inReal1[i]
	}
	return outReal
}

// Sub - Vector arithmetic subtraction
func Sub(inReal0 []float64, inReal1 []float64) []float64 {
	outReal := make([]float64, len(inReal0))
	for i := 0; i < len(inReal0); i++ {
		outReal[i] = inReal0[i] - inReal1[i]
	}
	return outReal
}

// Sum - Vector summation
func Sum(inReal []float64, optInTimePeriod int) []float64 {

	outReal := make([]float64, len(inReal))

	lookbackTotal := optInTimePeriod - 1
	startIdx := lookbackTotal
	periodTotal := 0.0
	trailingIdx := startIdx - lookbackTotal
	i := trailingIdx
	if optInTimePeriod > 1 {
		for i < startIdx {
			periodTotal += inReal[i]
			i++
		}
	}
	outIdx := startIdx
	for i < len(inReal) {
		periodTotal += inReal[i]
		tempReal := periodTotal
		periodTotal -= inReal[trailingIdx]
		outReal[outIdx] = tempReal
		i++
		trailingIdx++
		outIdx++
	}

	return outReal
}
