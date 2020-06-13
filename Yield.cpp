//============================================================================
//
// Author      :Mekkoaui Abdelhakim
// Version     : 1.0
// Description : yield method
//============================================================================

#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

class Bond {
public:
	Bond(double F, double issue_date, int num_periods, int freq,
			const vector<double> &c) {
		initValues(F, issue_date, num_periods, freq);
		setCoupons(c);
	}

	Bond(double F, double issue_date, int num_periods, int freq,
			const double c) {
		initValues(F, issue_date, num_periods, freq);
		setFlatCoupons(c);
	}

	~Bond() {
		_Face = 0;
		_cpnFreq = 0;
		_issue = 0;
		_maturity = 0;
		_cpnAmt.clear();
		_cpnDate.clear();
	}

	void setFlatCoupons(double c) {
		if (c < 0)
			c = 0.0;
		fill(_cpnAmt.begin(), _cpnAmt.end(), c);
	}

	void setCoupons(const vector<double> &c) {
		int i;
		for (i = 0; i < c.size(); i++) {
			if (c[i] >= 0.0)
				_cpnAmt[i] = c[i];
			else
				_cpnAmt[i] = 0.0;
		}
		if (c.size() < _numCpnPeriods) {
			for (int j = i; j < _cpnAmt.size(); j++) {
				if (c.back() >= 0.0)
					_cpnAmt[j] = c.back();
				else
					_cpnAmt[j] = 0.0;
			}
		}
	}

	double FairValue(double t0, double y) const {
		double B = 0;
		double dummy1 = 0; //Macaulay duration
		double dummy2 = 0; //Modified duration
		FV_duration(t0, y, B, dummy1, dummy2);
		return B;
	}

	double maturity() const {
		return _maturity;
	}
	double issue() const {
		return _issue;
	}

	int FV_duration(double t0, double y, double &B, double &Mac_dur,
			double &mod_dur) const {

		const double tol = 1.0e-6;
		double y_decimal = y * 0.01;
		if (t0 < issue() || t0 >= maturity())
			return 1;
		int i;
		for (i = 0; i < _numCpnPeriods - 1; i++) {

			//If we have t0 greater than the coupons dates because we do not have to pay coupons in the past.
			if (_cpnDate[i] < t0) {
				continue;
			}
			//without Face value
			if (_cpnDate[i] >= t0 + tol) {
				B = B
						+ (_cpnAmt[i] / _cpnFreq)
								/ pow((1 + (y_decimal / _cpnFreq)),
										_cpnFreq * (_cpnDate[i] - t0));
				Mac_dur = Mac_dur
						+ (_cpnDate[i] - t0) * (_cpnAmt[i] / _cpnFreq)
								/ pow((1 + (y_decimal / _cpnFreq)),
										_cpnFreq * (_cpnDate[i] - t0));
			}
		}
		//with Face value
		if (_cpnDate[i] >= t0 + tol) {
			B = B
					+ (_Face + (double(_cpnAmt[i]) / _cpnFreq))
							/ pow((1 + (y_decimal / _cpnFreq)),
									_cpnFreq * (_cpnDate[i] - t0));
			if (y_decimal > 0.0) {
				Mac_dur = (Mac_dur
						+ (_cpnDate[i] - t0)
								* (_Face + (double(_cpnAmt[i]) / _cpnFreq))
								/ pow((1 + (y_decimal / _cpnFreq)),
										_cpnFreq * (_cpnDate[i] - t0))) / B;
			} else
				Mac_dur = maturity() - t0;
		}
		mod_dur = Mac_dur / (1 + (y_decimal / _cpnFreq));

		return 0;
	}

private:
// data
	double _Face;
	double _issue;
	double _maturity;
	int _cpnFreq;
	int _numCpnPeriods;
	vector<double> _cpnAmt;
	vector<double> _cpnDate;

	void initValues(double F, double issue_date, int num_periods, int freq) {
		//Set Face = F. Impose the condition Face >= 0 in the constructor.
		if (F < 0)
			_Face = 0;
		else
			_Face = F;
		//Set the coupon frequency cpnFreq = freq.
		//Impose the condition cpnFreq >= 1 (at least one per year) in the constructor.
		if (freq >= 1)
			_cpnFreq = freq;
		else
			_cpnFreq = 1;
		/*Set the number of coupon periods numCpnPeriods = num periods.
		 Impose the condition numCpnPeriods >= 1 (at least one coupon period) in the
		 constructor.*/
		if (num_periods >= 1)
			_numCpnPeriods = num_periods;
		else
			_numCpnPeriods = 1;

		_issue = issue_date;

		_maturity = _issue + _numCpnPeriods / _cpnFreq;

		_cpnAmt.resize(_numCpnPeriods);
		_cpnDate.resize(_numCpnPeriods);

		//Calculate the coupon date.

		for (int i = 0; i < _numCpnPeriods; i++) {
			_cpnDate[i] = _issue + ((1.0 / _cpnFreq) * (i + 1));
		}
	}
};

int yield(double &y, int &num_iter, const Bond &bond, double B_target,
		double t0, double tol = 1.0e-4, int max_iter = 100) {

	y = 0;
	double y_low = 0.0;
	double B;

	if (B_target <= 0.0 || t0 < bond.issue() || t0 >= bond.maturity())
		return 1;

	double B_y_low = bond.FairValue(t0, y_low);
	double diff_B_y_low = B_y_low - B_target;

	if (fabs(diff_B_y_low) <= tol) {
		y = y_low;
		return 0;
	}

	double y_high = 100.0;
	double B_y_high = bond.FairValue(t0, y_high);
	double diff_B_y_high = B_y_high - B_target;

	if (fabs(diff_B_y_high) <= tol) {
		y = y_high;
		return 0;
	}

	if (diff_B_y_low * diff_B_y_high > 0) {
		y = 0;
		return 1;
	}

	double diff_B;

	for (num_iter = 1; num_iter < max_iter; ++num_iter) {
		y = (y_low + y_high) / 2;
		B = bond.FairValue(t0, y);
		diff_B = B - B_target;
		if (fabs(diff_B) <= tol)
			return 0;
		if (diff_B * diff_B_y_low > 0.0)
			y_low = y;
		else
			y_high = y;

		if (fabs(y_high - y_low) <= tol)
			return 0;

	}
	y = 0;
	return 0;
}

int main() {
	double b_duration = 0;
	double macaulay_duration = 0;
	double modified_duration = 0;
	int num_iter = 0;

	//Values to change
	//double c = 5; //Using a flat coupon
	//vector<double> c = {4.1,4.2,4.3,4.4}; //Using different coupons
	vector<double> c;
		c.push_back(2);
	  c.push_back(4);
	  c.push_back(6);
		c.push_back(8);
	//	c.push_back(5);
		//c.push_back(1);
	//	c.push_back(6);
		//c.push_back(6);//*/

	double face = 100.00;
	double issue = -0.1;
	int period = 4;
	int freq = 2;
	double t0 = 0.2;
	double B_target = 100;
	double y = 5.848;

	Bond b(face, issue, period, freq, c);

	//Results
	cout << "The value of FV= " << b.FairValue(t0, y) << endl;
	b.FV_duration(t0, y, b_duration, macaulay_duration, modified_duration);
	cout << "The Macaulay value is: " <<macaulay_duration << endl;
	cout << "The Modified duration value is: " <<modified_duration << endl <<endl;

	yield(y, num_iter, b, B_target, t0);
	cout << "Yield from bond price: " << y << endl;
	cout << "Number of iterations: " << num_iter;
}
