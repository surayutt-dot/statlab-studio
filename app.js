/**
 * StatLab Studio - Interactive Statistics Teaching Web Application
 * Features:
 * - Independent Variables (IV: Categorical grouping)
 * - Dependent Variables (DV: 5-Level Likert Rating Scale: 1=น้อยที่สุด to 5=มากที่สุด)
 * - Pure HTML Math Symbols (No raw LaTeX dollar strings)
 * - Frequency, Percentage, Mean, S.D., 5-Level Interpretation, One-Sample t-test,
 *   Independent Samples t-test (with Levene's test), One-Way ANOVA (with Post-hoc LSD),
 *   and Practice Table-Reading Quizzes.
 */

// ============================================================================
// 1. 5-LEVEL LIKERT SCALE DEFINITIONS & INTERPRETATION HELPER
// ============================================================================

const LIKERT_LEVELS = [
  { val: 5, label: "5 - มากที่สุด", short: "มากที่สุด", badgeClass: "bg-emerald-100 text-emerald-800 border-emerald-200" },
  { val: 4, label: "4 - มาก", short: "มาก", badgeClass: "bg-teal-100 text-teal-800 border-teal-200" },
  { val: 3, label: "3 - ปานกลาง", short: "ปานกลาง", badgeClass: "bg-amber-100 text-amber-800 border-amber-200" },
  { val: 2, label: "2 - น้อย", short: "น้อย", badgeClass: "bg-orange-100 text-orange-800 border-orange-200" },
  { val: 1, label: "1 - น้อยที่สุด", short: "น้อยที่สุด", badgeClass: "bg-rose-100 text-rose-800 border-rose-200" }
];

function interpretMean5Level(meanVal) {
  if (isNaN(meanVal) || meanVal <= 0) {
    return { text: "-", badge: `<span class="text-slate-400">-</span>` };
  }
  if (meanVal >= 4.21) {
    return {
      text: "มากที่สุด",
      badge: `<span class="inline-block bg-emerald-100 text-emerald-800 border border-emerald-200 px-2 py-0.5 rounded text-xs font-semibold">มากที่สุด</span>`
    };
  }
  if (meanVal >= 3.41) {
    return {
      text: "มาก",
      badge: `<span class="inline-block bg-teal-100 text-teal-800 border border-teal-200 px-2 py-0.5 rounded text-xs font-semibold">มาก</span>`
    };
  }
  if (meanVal >= 2.61) {
    return {
      text: "ปานกลาง",
      badge: `<span class="inline-block bg-amber-100 text-amber-800 border border-amber-200 px-2 py-0.5 rounded text-xs font-semibold">ปานกลาง</span>`
    };
  }
  if (meanVal >= 1.81) {
    return {
      text: "น้อย",
      badge: `<span class="inline-block bg-orange-100 text-orange-800 border border-orange-200 px-2 py-0.5 rounded text-xs font-semibold">น้อย</span>`
    };
  }
  return {
    text: "น้อยที่สุด",
    badge: `<span class="inline-block bg-rose-100 text-rose-800 border border-rose-200 px-2 py-0.5 rounded text-xs font-semibold">น้อยที่สุด</span>`
  };
}

// ============================================================================
// 2. MATHEMATICAL & STATISTICAL ENGINE (Exact p-values via Incomplete Beta)
// ============================================================================

const StatMath = {
  mean(arr) {
    if (!arr || arr.length === 0) return 0;
    return arr.reduce((a, b) => a + b, 0) / arr.length;
  },

  variance(arr) {
    if (!arr || arr.length < 2) return 0;
    const m = this.mean(arr);
    return arr.reduce((acc, val) => acc + Math.pow(val - m, 2), 0) / (arr.length - 1);
  },

  sd(arr) {
    return Math.sqrt(this.variance(arr));
  },

  se(arr) {
    if (!arr || arr.length === 0) return 0;
    return this.sd(arr) / Math.sqrt(arr.length);
  },

  logGamma(z) {
    const c = [
      76.18009172947146,
      -86.50532032941677,
      24.01409824083091,
      -1.231739572450155,
      0.001208650973866179,
      -0.000005395239384953
    ];
    let y = z;
    let x = z;
    let tmp = x + 5.5;
    tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) {
      y += 1;
      ser += c[j] / y;
    }
    return -tmp + Math.log((2.5066282746310005 * ser) / x);
  },

  betacf(a, b, x) {
    const MAXIT = 200;
    const EPS = 3.0e-9;
    const FPMIN = 1.0e-30;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let c = 1.0;
    let d = 1.0 - (qab * x) / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    let h = d;
    for (let m = 1; m <= MAXIT; m++) {
      let m2 = 2 * m;
      let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
      d = 1.0 + aa * d;
      if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1.0 + aa / c;
      if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1.0 / d;
      h *= d * c;
      aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
      d = 1.0 + aa * d;
      if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1.0 + aa / c;
      if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1.0 / d;
      let del = d * c;
      h *= del;
      if (Math.abs(del - 1.0) < EPS) break;
    }
    return h;
  },

  betai(a, b, x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const bt = Math.exp(
      this.logGamma(a + b) -
        this.logGamma(a) -
        this.logGamma(b) +
        a * Math.log(x) +
        b * Math.log(1.0 - x)
    );
    if (x < (a + 1.0) / (a + b + 2.0)) {
      return (bt * this.betacf(a, b, x)) / a;
    } else {
      return 1.0 - (bt * this.betacf(b, a, 1.0 - x)) / b;
    }
  },

  tProbTwoTailed(tVal, df) {
    if (!isFinite(tVal) || df <= 0) return 1.0;
    const x = df / (df + tVal * tVal);
    return Math.min(1, Math.max(0, this.betai(0.5 * df, 0.5, x)));
  },

  tCritical(alpha, df) {
    if (df <= 0) return 1.96;
    let low = 0;
    let high = 50;
    for (let i = 0; i < 60; i++) {
      let mid = (low + high) / 2;
      let p = this.tProbTwoTailed(mid, df);
      if (p > alpha) low = mid;
      else high = mid;
    }
    return (low + high) / 2;
  },

  fProbRightTailed(fVal, df1, df2) {
    if (!isFinite(fVal) || fVal <= 0 || df1 <= 0 || df2 <= 0) return 1.0;
    const x = df2 / (df2 + df1 * fVal);
    return Math.min(1, Math.max(0, this.betai(0.5 * df2, 0.5 * df1, x)));
  },

  oneSampleTTest(values, mu0, alpha = 0.05) {
    const n = values.length;
    const mean = this.mean(values);
    const sd = this.sd(values);
    const se = sd / Math.sqrt(n);
    const diff = mean - mu0;
    const df = n - 1;
    const t = se > 0 ? diff / se : 0;
    const p = this.tProbTwoTailed(t, df);
    const tCrit = this.tCritical(alpha, df);
    const ciLow = diff - tCrit * se;
    const ciHigh = diff + tCrit * se;
    return { n, mean, sd, se, mu0, diff, t, df, p, ciLow, ciHigh };
  },

  independentTTest(group1Vals, group2Vals) {
    const n1 = group1Vals.length;
    const n2 = group2Vals.length;
    const m1 = this.mean(group1Vals);
    const m2 = this.mean(group2Vals);
    const v1 = this.variance(group1Vals);
    const v2 = this.variance(group2Vals);
    const sd1 = Math.sqrt(v1);
    const sd2 = Math.sqrt(v2);
    const se1 = sd1 / Math.sqrt(n1);
    const se2 = sd2 / Math.sqrt(n2);
    const meanDiff = m1 - m2;

    const z1 = group1Vals.map((y) => Math.abs(y - m1));
    const z2 = group2Vals.map((y) => Math.abs(y - m2));
    const zm1 = this.mean(z1);
    const zm2 = this.mean(z2);
    const zAll = [...z1, ...z2];
    const zmTotal = this.mean(zAll);
    const N = n1 + n2;
    const ssBetweenZ = n1 * Math.pow(zm1 - zmTotal, 2) + n2 * Math.pow(zm2 - zmTotal, 2);
    const ssWithinZ =
      z1.reduce((acc, val) => acc + Math.pow(val - zm1, 2), 0) +
      z2.reduce((acc, val) => acc + Math.pow(val - zm2, 2), 0);
    const dfLevene2 = N - 2;
    const msWithinZ = dfLevene2 > 0 ? ssWithinZ / dfLevene2 : 0;
    const leveneF = msWithinZ > 0 ? ssBetweenZ / msWithinZ : 0;
    const leveneP = this.fProbRightTailed(leveneF, 1, dfLevene2);

    const dfEqual = n1 + n2 - 2;
    const sp2 = dfEqual > 0 ? ((n1 - 1) * v1 + (n2 - 1) * v2) / dfEqual : 0;
    const seDiffEqual = Math.sqrt(sp2 * (1 / n1 + 1 / n2));
    const tEqual = seDiffEqual > 0 ? meanDiff / seDiffEqual : 0;
    const pEqual = this.tProbTwoTailed(tEqual, dfEqual);

    const seDiffUnequal = Math.sqrt(v1 / n1 + v2 / n2);
    const tUnequal = seDiffUnequal > 0 ? meanDiff / seDiffUnequal : 0;
    const numWelch = Math.pow(v1 / n1 + v2 / n2, 2);
    const denWelch =
      Math.pow(v1 / n1, 2) / Math.max(1, n1 - 1) +
      Math.pow(v2 / n2, 2) / Math.max(1, n2 - 1);
    const dfUnequal = denWelch > 0 ? numWelch / denWelch : dfEqual;
    const pUnequal = this.tProbTwoTailed(tUnequal, dfUnequal);

    const equalAssumed = leveneP > 0.05;

    return {
      g1: { n: n1, mean: m1, sd: sd1, se: se1 },
      g2: { n: n2, mean: m2, sd: sd2, se: se2 },
      meanDiff,
      leveneF,
      leveneP,
      equalAssumed,
      equalRow: { t: tEqual, df: dfEqual, p: pEqual, seDiff: seDiffEqual },
      unequalRow: { t: tUnequal, df: dfUnequal, p: pUnequal, seDiff: seDiffUnequal }
    };
  },

  oneWayAnova(groupsObj) {
    const validGroups = groupsObj.filter((g) => g.values.length > 0);
    const k = validGroups.length;
    const allValues = validGroups.flatMap((g) => g.values);
    const N = allValues.length;
    const grandMean = this.mean(allValues);

    const descriptives = validGroups.map((g) => {
      const n = g.values.length;
      const mean = this.mean(g.values);
      const sd = this.sd(g.values);
      const se = this.se(g.values);
      const min = Math.min(...g.values);
      const max = Math.max(...g.values);
      return { name: g.name, n, mean, sd, se, min, max, values: g.values };
    });

    let ssBetween = 0;
    let ssWithin = 0;
    descriptives.forEach((d) => {
      ssBetween += d.n * Math.pow(d.mean - grandMean, 2);
      d.values.forEach((v) => {
        ssWithin += Math.pow(v - d.mean, 2);
      });
    });
    const ssTotal = ssBetween + ssWithin;

    const dfBetween = Math.max(1, k - 1);
    const dfWithin = Math.max(1, N - k);
    const dfTotal = Math.max(1, N - 1);

    const msBetween = ssBetween / dfBetween;
    const msWithin = ssWithin / dfWithin;
    const fVal = msWithin > 0 ? msBetween / msWithin : 0;
    const pVal = this.fProbRightTailed(fVal, dfBetween, dfWithin);

    const posthoc = [];
    for (let i = 0; i < descriptives.length; i++) {
      for (let j = i + 1; j < descriptives.length; j++) {
        const gI = descriptives[i];
        const gJ = descriptives[j];
        const diff = gI.mean - gJ.mean;
        const sePair = Math.sqrt(msWithin * (1 / gI.n + 1 / gJ.n));
        const tPair = sePair > 0 ? diff / sePair : 0;
        const pPair = this.tProbTwoTailed(tPair, dfWithin);
        posthoc.push({
          groupI: gI.name,
          groupJ: gJ.name,
          meanI: gI.mean,
          meanJ: gJ.mean,
          diff,
          se: sePair,
          t: tPair,
          p: pPair
        });
      }
    }

    return {
      N,
      k,
      grandMean,
      grandSd: this.sd(allValues),
      descriptives,
      ssBetween,
      ssWithin,
      ssTotal,
      dfBetween,
      dfWithin,
      dfTotal,
      msBetween,
      msWithin,
      fVal,
      pVal,
      posthoc
    };
  }
};

// ============================================================================
// 3. PRESET TEACHING DATASETS (5-Level Rating Scale DVs: 1 to 5)
// ============================================================================

const PRESETS = {
  student_cafe: {
    title: "ความพึงพอใจของนักศึกษาต่อสิ่งอำนวยความสะดวกและบริการของมหาวิทยาลัย",
    description:
      "ศึกษาเปรียบเทียบระดับความพึงพอใจ 5 ระดับ (1=น้อยที่สุด ถึง 5=มากที่สุด) จำแนกตามเพศ ลักษณะที่พักอาศัย และชั้นปี",
    ivs: [
      { id: "gender", label: "เพศสภาพ", options: ["ชาย", "หญิง"] },
      { id: "residence", label: "ลักษณะที่พัก", options: ["หอพักใน", "หอพักนอก", "บ้านตนเอง"] },
      { id: "year", label: "ชั้นปี", options: ["ปี 1", "ปี 2", "ปี 3", "ปี 4"] }
    ],
    dvs: [
      { id: "canteen_sat", label: "ความพึงพอใจด้านคุณภาพโรงอาหาร", unit: "ระดับ (1-5)", defaultTestVal: 3.5 },
      { id: "wifi_sat", label: "ความพึงพอใจด้านสัญญาณอินเทอร์เน็ต Wi-Fi", unit: "ระดับ (1-5)", defaultTestVal: 3.5 },
      { id: "library_sat", label: "ความพึงพอใจด้านพื้นที่อ่านหนังสือห้องสมุด", unit: "ระดับ (1-5)", defaultTestVal: 3.5 }
    ],
    data: [
      { gender: "หญิง", residence: "หอพักนอก", year: "ปี 2", canteen_sat: 4, wifi_sat: 3, library_sat: 5 },
      { gender: "ชาย", residence: "หอพักใน", year: "ปี 1", canteen_sat: 3, wifi_sat: 2, library_sat: 4 },
      { gender: "หญิง", residence: "หอพักนอก", year: "ปี 3", canteen_sat: 5, wifi_sat: 4, library_sat: 5 },
      { gender: "ชาย", residence: "บ้านตนเอง", year: "ปี 2", canteen_sat: 3, wifi_sat: 3, library_sat: 3 },
      { gender: "หญิง", residence: "หอพักใน", year: "ปี 1", canteen_sat: 4, wifi_sat: 2, library_sat: 4 },
      { gender: "ชาย", residence: "หอพักนอก", year: "ปี 4", canteen_sat: 4, wifi_sat: 4, library_sat: 4 },
      { gender: "หญิง", residence: "บ้านตนเอง", year: "ปี 3", canteen_sat: 4, wifi_sat: 3, library_sat: 5 },
      { gender: "ชาย", residence: "หอพักใน", year: "ปี 1", canteen_sat: 2, wifi_sat: 2, library_sat: 3 },
      { gender: "หญิง", residence: "หอพักนอก", year: "ปี 4", canteen_sat: 5, wifi_sat: 4, library_sat: 5 },
      { gender: "ชาย", residence: "หอพักนอก", year: "ปี 2", canteen_sat: 4, wifi_sat: 3, library_sat: 4 },
      { gender: "หญิง", residence: "หอพักใน", year: "ปี 2", canteen_sat: 4, wifi_sat: 2, library_sat: 4 },
      { gender: "ชาย", residence: "บ้านตนเอง", year: "ปี 3", canteen_sat: 3, wifi_sat: 3, library_sat: 3 },
      { gender: "หญิง", residence: "หอพักนอก", year: "ปี 3", canteen_sat: 5, wifi_sat: 4, library_sat: 5 },
      { gender: "ชาย", residence: "หอพักใน", year: "ปี 1", canteen_sat: 2, wifi_sat: 1, library_sat: 3 },
      { gender: "หญิง", residence: "บ้านตนเอง", year: "ปี 4", canteen_sat: 4, wifi_sat: 3, library_sat: 4 },
      { gender: "ชาย", residence: "หอพักนอก", year: "ปี 3", canteen_sat: 4, wifi_sat: 4, library_sat: 4 },
      { gender: "หญิง", residence: "หอพักใน", year: "ปี 1", canteen_sat: 3, wifi_sat: 2, library_sat: 4 },
      { gender: "ชาย", residence: "บ้านตนเอง", year: "ปี 4", canteen_sat: 3, wifi_sat: 3, library_sat: 3 },
      { gender: "หญิง", residence: "หอพักนอก", year: "ปี 2", canteen_sat: 5, wifi_sat: 4, library_sat: 5 },
      { gender: "ชาย", residence: "หอพักใน", year: "ปี 2", canteen_sat: 3, wifi_sat: 2, library_sat: 3 },
      { gender: "หญิง", residence: "หอพักนอก", year: "ปี 4", canteen_sat: 5, wifi_sat: 5, library_sat: 5 },
      { gender: "ชาย", residence: "หอพักนอก", year: "ปี 1", canteen_sat: 4, wifi_sat: 3, library_sat: 4 },
      { gender: "หญิง", residence: "บ้านตนเอง", year: "ปี 2", canteen_sat: 4, wifi_sat: 3, library_sat: 4 },
      { gender: "ชาย", residence: "หอพักใน", year: "ปี 3", canteen_sat: 3, wifi_sat: 2, library_sat: 4 },
      { gender: "หญิง", residence: "หอพักใน", year: "ปี 3", canteen_sat: 4, wifi_sat: 3, library_sat: 4 },
      { gender: "ชาย", residence: "บ้านตนเอง", year: "ปี 1", canteen_sat: 2, wifi_sat: 3, library_sat: 3 },
      { gender: "หญิง", residence: "หอพักนอก", year: "ปี 1", canteen_sat: 4, wifi_sat: 4, library_sat: 5 },
      { gender: "ชาย", residence: "หอพักนอก", year: "ปี 4", canteen_sat: 4, wifi_sat: 4, library_sat: 4 },
      { gender: "หญิง", residence: "บ้านตนเอง", year: "ปี 3", canteen_sat: 4, wifi_sat: 3, library_sat: 5 },
      { gender: "ชาย", residence: "หอพักใน", year: "ปี 2", canteen_sat: 3, wifi_sat: 2, library_sat: 3 },
      { gender: "หญิง", residence: "หอพักนอก", year: "ปี 4", canteen_sat: 5, wifi_sat: 4, library_sat: 5 },
      { gender: "ชาย", residence: "บ้านตนเอง", year: "ปี 2", canteen_sat: 3, wifi_sat: 3, library_sat: 3 },
      { gender: "หญิง", residence: "หอพักใน", year: "ปี 2", canteen_sat: 4, wifi_sat: 2, library_sat: 4 },
      { gender: "ชาย", residence: "หอพักนอก", year: "ปี 3", canteen_sat: 4, wifi_sat: 4, library_sat: 4 },
      { gender: "หญิง", residence: "บ้านตนเอง", year: "ปี 1", canteen_sat: 3, wifi_sat: 3, library_sat: 4 },
      { gender: "ชาย", residence: "หอพักใน", year: "ปี 4", canteen_sat: 3, wifi_sat: 2, library_sat: 4 }
    ]
  },

  teaching_method: {
    title: "การประเมินผลการจัดการเรียนรู้วิชาสถิติ 3 รูปแบบ (มาตรวัด 5 ระดับ)",
    description:
      "เปรียบเทียบวิธีสอน 3 รูปแบบ และรูปแบบการทำกิจกรรม ต่อระดับความเข้าใจเนื้อหาและระดับความพึงพอใจในการเรียน (1-5 คะแนน)",
    ivs: [
      { id: "method", label: "รูปแบบวิธีการสอน", options: ["บรรยายปกติ", "ใช้เกม (Gamification)", "ปฏิบัติการ (Active Lab)"] },
      { id: "work_mode", label: "ลักษณะการทำแบบฝึกหัด", options: ["ทำเดี่ยว", "ทำกลุ่มคู่"] }
    ],
    dvs: [
      { id: "understand_lvl", label: "ระดับความเข้าใจในเนื้อหาสถิติ", unit: "ระดับ (1-5)", defaultTestVal: 3.5 },
      { id: "enjoy_lvl", label: "ระดับความสนุกและพึงพอใจต่อการเรียน", unit: "ระดับ (1-5)", defaultTestVal: 3.5 }
    ],
    data: [
      { method: "บรรยายปกติ", work_mode: "ทำเดี่ยว", understand_lvl: 3, enjoy_lvl: 2 },
      { method: "บรรยายปกติ", work_mode: "ทำเดี่ยว", understand_lvl: 3, enjoy_lvl: 3 },
      { method: "บรรยายปกติ", work_mode: "ทำกลุ่มคู่", understand_lvl: 3, enjoy_lvl: 3 },
      { method: "บรรยายปกติ", work_mode: "ทำเดี่ยว", understand_lvl: 2, enjoy_lvl: 2 },
      { method: "บรรยายปกติ", work_mode: "ทำกลุ่มคู่", understand_lvl: 3, enjoy_lvl: 3 },
      { method: "บรรยายปกติ", work_mode: "ทำเดี่ยว", understand_lvl: 3, enjoy_lvl: 2 },
      { method: "บรรยายปกติ", work_mode: "ทำกลุ่มคู่", understand_lvl: 4, enjoy_lvl: 3 },
      { method: "บรรยายปกติ", work_mode: "ทำเดี่ยว", understand_lvl: 2, enjoy_lvl: 2 },
      { method: "บรรยายปกติ", work_mode: "ทำกลุ่มคู่", understand_lvl: 3, enjoy_lvl: 3 },
      { method: "บรรยายปกติ", work_mode: "ทำเดี่ยว", understand_lvl: 3, enjoy_lvl: 2 },
      { method: "ใช้เกม (Gamification)", work_mode: "ทำกลุ่มคู่", understand_lvl: 4, enjoy_lvl: 5 },
      { method: "ใช้เกม (Gamification)", work_mode: "ทำเดี่ยว", understand_lvl: 4, enjoy_lvl: 4 },
      { method: "ใช้เกม (Gamification)", work_mode: "ทำกลุ่มคู่", understand_lvl: 4, enjoy_lvl: 5 },
      { method: "ใช้เกม (Gamification)", work_mode: "ทำเดี่ยว", understand_lvl: 3, enjoy_lvl: 4 },
      { method: "ใช้เกม (Gamification)", work_mode: "ทำกลุ่มคู่", understand_lvl: 4, enjoy_lvl: 5 },
      { method: "ใช้เกม (Gamification)", work_mode: "ทำเดี่ยว", understand_lvl: 4, enjoy_lvl: 4 },
      { method: "ใช้เกม (Gamification)", work_mode: "ทำกลุ่มคู่", understand_lvl: 4, enjoy_lvl: 5 },
      { method: "ใช้เกม (Gamification)", work_mode: "ทำเดี่ยว", understand_lvl: 4, enjoy_lvl: 4 },
      { method: "ใช้เกม (Gamification)", work_mode: "ทำกลุ่มคู่", understand_lvl: 5, enjoy_lvl: 5 },
      { method: "ใช้เกม (Gamification)", work_mode: "ทำเดี่ยว", understand_lvl: 4, enjoy_lvl: 4 },
      { method: "ปฏิบัติการ (Active Lab)", work_mode: "ทำกลุ่มคู่", understand_lvl: 5, enjoy_lvl: 5 },
      { method: "ปฏิบัติการ (Active Lab)", work_mode: "ทำเดี่ยว", understand_lvl: 4, enjoy_lvl: 4 },
      { method: "ปฏิบัติการ (Active Lab)", work_mode: "ทำกลุ่มคู่", understand_lvl: 5, enjoy_lvl: 5 },
      { method: "ปฏิบัติการ (Active Lab)", work_mode: "ทำเดี่ยว", understand_lvl: 5, enjoy_lvl: 4 },
      { method: "ปฏิบัติการ (Active Lab)", work_mode: "ทำกลุ่มคู่", understand_lvl: 5, enjoy_lvl: 5 },
      { method: "ปฏิบัติการ (Active Lab)", work_mode: "ทำเดี่ยว", understand_lvl: 4, enjoy_lvl: 4 },
      { method: "ปฏิบัติการ (Active Lab)", work_mode: "ทำกลุ่มคู่", understand_lvl: 5, enjoy_lvl: 5 },
      { method: "ปฏิบัติการ (Active Lab)", work_mode: "ทำเดี่ยว", understand_lvl: 5, enjoy_lvl: 4 },
      { method: "ปฏิบัติการ (Active Lab)", work_mode: "ทำกลุ่มคู่", understand_lvl: 5, enjoy_lvl: 5 },
      { method: "ปฏิบัติการ (Active Lab)", work_mode: "ทำกลุ่มคู่", understand_lvl: 5, enjoy_lvl: 5 }
    ]
  },

  health_sleep: {
    title: "การสำรวจพฤติกรรมสุขภาพ คุณภาพการนอนหลับ และความสดชื่น (5 ระดับ)",
    description:
      "สำรวจความถี่ในการออกกำลังกายและการดื่มกาแฟ ที่มีผลต่อระดับคุณภาพการนอนหลับและระดับความสดชื่นระหว่างวัน (1-5 คะแนน)",
    ivs: [
      { id: "coffee", label: "การดื่มกาแฟเป็นประจำ", options: ["ดื่มทุกวัน", "ไม่ดื่ม"] },
      { id: "exercise", label: "ความถี่การออกกำลังกาย", options: ["ไม่ออกกำลังกาย", "1-2 ครั้ง/สัปดาห์", "3 ครั้งขึ้นไป/สัปดาห์"] }
    ],
    dvs: [
      { id: "sleep_quality", label: "ระดับคุณภาพการนอนหลับสนิท", unit: "ระดับ (1-5)", defaultTestVal: 3.5 },
      { id: "freshness", label: "ระดับความสดชื่นกระปรี้กระเปร่า", unit: "ระดับ (1-5)", defaultTestVal: 3.5 }
    ],
    data: [
      { coffee: "ดื่มทุกวัน", exercise: "ไม่ออกกำลังกาย", sleep_quality: 2, freshness: 2 },
      { coffee: "ดื่มทุกวัน", exercise: "ไม่ออกกำลังกาย", sleep_quality: 2, freshness: 3 },
      { coffee: "ไม่ดื่ม", exercise: "1-2 ครั้ง/สัปดาห์", sleep_quality: 4, freshness: 4 },
      { coffee: "ไม่ดื่ม", exercise: "3 ครั้งขึ้นไป/สัปดาห์", sleep_quality: 5, freshness: 5 },
      { coffee: "ดื่มทุกวัน", exercise: "1-2 ครั้ง/สัปดาห์", sleep_quality: 3, freshness: 3 },
      { coffee: "ไม่ดื่ม", exercise: "ไม่ออกกำลังกาย", sleep_quality: 3, freshness: 2 },
      { coffee: "ดื่มทุกวัน", exercise: "3 ครั้งขึ้นไป/สัปดาห์", sleep_quality: 4, freshness: 4 },
      { coffee: "ไม่ดื่ม", exercise: "3 ครั้งขึ้นไป/สัปดาห์", sleep_quality: 5, freshness: 5 },
      { coffee: "ดื่มทุกวัน", exercise: "ไม่ออกกำลังกาย", sleep_quality: 1, freshness: 2 },
      { coffee: "ไม่ดื่ม", exercise: "1-2 ครั้ง/สัปดาห์", sleep_quality: 4, freshness: 4 },
      { coffee: "ดื่มทุกวัน", exercise: "1-2 ครั้ง/สัปดาห์", sleep_quality: 3, freshness: 3 },
      { coffee: "ไม่ดื่ม", exercise: "3 ครั้งขึ้นไป/สัปดาห์", sleep_quality: 5, freshness: 5 },
      { coffee: "ดื่มทุกวัน", exercise: "ไม่ออกกำลังกาย", sleep_quality: 2, freshness: 2 },
      { coffee: "ไม่ดื่ม", exercise: "ไม่ออกกำลังกาย", sleep_quality: 3, freshness: 3 },
      { coffee: "ดื่มทุกวัน", exercise: "3 ครั้งขึ้นไป/สัปดาห์", sleep_quality: 4, freshness: 4 },
      { coffee: "ไม่ดื่ม", exercise: "1-2 ครั้ง/สัปดาห์", sleep_quality: 4, freshness: 4 },
      { coffee: "ดื่มทุกวัน", exercise: "ไม่ออกกำลังกาย", sleep_quality: 2, freshness: 2 },
      { coffee: "ไม่ดื่ม", exercise: "3 ครั้งขึ้นไป/สัปดาห์", sleep_quality: 5, freshness: 5 },
      { coffee: "ดื่มทุกวัน", exercise: "1-2 ครั้ง/สัปดาห์", sleep_quality: 3, freshness: 4 },
      { coffee: "ไม่ดื่ม", exercise: "ไม่ออกกำลังกาย", sleep_quality: 3, freshness: 2 },
      { coffee: "ดื่มทุกวัน", exercise: "3 ครั้งขึ้นไป/สัปดาห์", sleep_quality: 4, freshness: 4 },
      { coffee: "ไม่ดื่ม", exercise: "1-2 ครั้ง/สัปดาห์", sleep_quality: 4, freshness: 4 },
      { coffee: "ดื่มทุกวัน", exercise: "ไม่ออกกำลังกาย", sleep_quality: 2, freshness: 2 },
      { coffee: "ไม่ดื่ม", exercise: "3 ครั้งขึ้นไป/สัปดาห์", sleep_quality: 5, freshness: 5 },
      { coffee: "ดื่มทุกวัน", exercise: "1-2 ครั้ง/สัปดาห์", sleep_quality: 3, freshness: 3 },
      { coffee: "ไม่ดื่ม", exercise: "1-2 ครั้ง/สัปดาห์", sleep_quality: 4, freshness: 4 },
      { coffee: "ดื่มทุกวัน", exercise: "3 ครั้งขึ้นไป/สัปดาห์", sleep_quality: 4, freshness: 4 },
      { coffee: "ไม่ดื่ม", exercise: "ไม่ออกกำลังกาย", sleep_quality: 3, freshness: 3 },
      { coffee: "ดื่มทุกวัน", exercise: "1-2 ครั้ง/สัปดาห์", sleep_quality: 3, freshness: 3 },
      { coffee: "ไม่ดื่ม", exercise: "3 ครั้งขึ้นไป/สัปดาห์", sleep_quality: 5, freshness: 5 },
      { coffee: "ดื่มทุกวัน", exercise: "ไม่ออกกำลังกาย", sleep_quality: 2, freshness: 2 },
      { coffee: "ไม่ดื่ม", exercise: "1-2 ครั้ง/สัปดาห์", sleep_quality: 4, freshness: 4 }
    ]
  },

  custom_blank: {
    title: "หัวข้อการศึกษาของฉัน (คลิกปุ่มแก้ไขชื่อหัวข้อวิจัย)",
    description:
      "ให้นักศึกษากำหนดตัวแปรต้น (IV) และตัวแปรตาม 5 ระดับ (DV) ของกลุ่มตนเอง แล้วนำไปเก็บข้อมูลจากเพื่อนหรือกลุ่มตัวอย่างจริง",
    ivs: [
      { id: "gender", label: "เพศ", options: ["ชาย", "หญิง"] },
      { id: "group_var", label: "กลุ่มตัวอย่าง (เช่น ชั้นปี/สาขา)", options: ["กลุ่มที่ 1", "กลุ่มที่ 2", "กลุ่มที่ 3"] }
    ],
    dvs: [
      { id: "q1", label: "ข้อคำถามที่ 1 (ความพึงพอใจ/ความคิดเห็นด้านที่ 1)", unit: "ระดับ (1-5)", defaultTestVal: 3.5 },
      { id: "q2", label: "ข้อคำถามที่ 2 (ความพึงพอใจ/ความคิดเห็นด้านที่ 2)", unit: "ระดับ (1-5)", defaultTestVal: 3.5 }
    ],
    data: []
  }
};

// ============================================================================
// 4. APPLICATION STATE, LOCALSTORAGE AUTO-SAVE & HELPERS
// ============================================================================

const STORAGE_KEY = "statlab_studio_saved_state_v1";

let appState = {
  presetKey: "student_cafe",
  title: "",
  description: "",
  ivs: [], // [{ id, label, options: [] }]
  dvs: [], // [{ id, label, unit, defaultTestVal }]
  data: [], // Array of record objects with DV values 1..5
  currentQuiz: null
};

function saveStateToLocal() {
  try {
    const payload = {
      presetKey: appState.presetKey,
      title: appState.title,
      description: appState.description,
      ivs: appState.ivs,
      dvs: appState.dvs,
      data: appState.data
    };
    localStorage.setItem(STORAGE_KEY, JSON.stringify(payload));
    const badge = document.getElementById("autoSaveStatus");
    if (badge) {
      badge.innerHTML = `<span>💾</span> บันทึกอัตโนมัติแล้ว`;
    }
  } catch (err) {
    console.warn("Could not save to localStorage:", err);
  }
}

function loadStateFromLocal() {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return false;
    const parsed = JSON.parse(raw);
    if (!parsed || !Array.isArray(parsed.ivs) || !Array.isArray(parsed.dvs) || !Array.isArray(parsed.data)) {
      return false;
    }
    appState.presetKey = parsed.presetKey || "student_cafe";
    appState.title = parsed.title || "";
    appState.description = parsed.description || "";
    appState.ivs = parsed.ivs;
    appState.dvs = parsed.dvs;
    appState.data = parsed.data;
    const sel = document.getElementById("presetSelect");
    if (sel && PRESETS[appState.presetKey]) {
      sel.value = appState.presetKey;
    }
    return true;
  } catch (err) {
    return false;
  }
}

function fmt(num, digits = 2) {
  if (num === null || num === undefined || isNaN(num)) return "-";
  return Number(num).toFixed(digits);
}

function fmtP(p) {
  if (p === null || p === undefined || isNaN(p)) return "-";
  if (p < 0.001) return "< .001*";
  return p.toFixed(3) + (p < 0.05 ? "*" : "");
}

function fmtPRaw(p) {
  if (p === null || p === undefined || isNaN(p)) return "-";
  if (p < 0.001) return ".000";
  return p.toFixed(3);
}

function clampLikert(val) {
  const n = Math.round(Number(val));
  if (isNaN(n)) return 3;
  return Math.min(5, Math.max(1, n));
}

function loadPreset(presetKey) {
  const preset = PRESETS[presetKey] || PRESETS.student_cafe;
  appState.presetKey = presetKey;
  appState.title = preset.title;
  appState.description = preset.description;
  appState.ivs = JSON.parse(JSON.stringify(preset.ivs));
  appState.dvs = JSON.parse(JSON.stringify(preset.dvs));
  appState.data = JSON.parse(JSON.stringify(preset.data));

  const sel = document.getElementById("presetSelect");
  if (sel) sel.value = presetKey;

  saveStateToLocal();
  renderAll();
}

// ============================================================================
// 5. UI RENDERING: TAB 1 (VARIABLE DESIGN) & TAB 2 (DATA COLLECTION)
// ============================================================================

function renderAll() {
  renderProjectHeader();
  renderVariableCards();
  renderSurveyForm();
  renderRawDataTable();
  populateAnalysisSelectors();
  runAllAnalyses();
  generateQuiz();
}

function renderProjectHeader() {
  document.getElementById("projectTitleDisplay").textContent = `หัวข้อศึกษา: ${appState.title}`;
  document.getElementById("projectDescDisplay").textContent = appState.description;
  document.getElementById("badgeCount").textContent = appState.data.length;
  document.getElementById("totalRowsLabel").textContent = appState.data.length;
  document.getElementById("nextRespondentNum").textContent = appState.data.length + 1;
}

function renderVariableCards() {
  // Independent Variables (IVs)
  const ivContainer = document.getElementById("ivListContainer");
  ivContainer.innerHTML = "";
  appState.ivs.forEach((iv, idx) => {
    const groupCountBadge =
      iv.options.length === 2
        ? `<span class="text-[11px] bg-indigo-100 text-indigo-800 px-2 py-0.5 rounded font-medium">2 กลุ่ม (เหมาะกับ Independent t-test)</span>`
        : `<span class="text-[11px] bg-purple-100 text-purple-800 px-2 py-0.5 rounded font-medium">${iv.options.length} กลุ่ม (เหมาะกับ One-Way ANOVA)</span>`;

    const div = document.createElement("div");
    div.className = "bg-white p-3 rounded-lg border border-sky-200 shadow-sm flex items-start justify-between gap-2";
    div.innerHTML = `
      <div>
        <div class="flex flex-wrap items-center gap-2">
          <span class="font-heading font-semibold text-sm text-slate-800">${idx + 1}. ${iv.label}</span>
          ${groupCountBadge}
        </div>
        <div class="flex flex-wrap gap-1.5 mt-1.5">
          ${iv.options
            .map(
              (opt) =>
                `<span class="text-xs bg-slate-100 text-slate-700 px-2 py-0.5 rounded border border-slate-200">${opt}</span>`
            )
            .join("")}
        </div>
      </div>
      ${
        appState.ivs.length > 1
          ? `<button data-remove-iv="${iv.id}" class="text-xs text-rose-500 hover:text-rose-700 px-1.5 py-0.5" title="ลบตัวแปร">✕</button>`
          : ""
      }
    `;
    ivContainer.appendChild(div);
  });

  // Dependent Variables (DVs - 5 Levels)
  const dvContainer = document.getElementById("dvListContainer");
  dvContainer.innerHTML = "";
  appState.dvs.forEach((dv, idx) => {
    const div = document.createElement("div");
    div.className = "bg-white p-3 rounded-lg border border-emerald-200 shadow-sm flex items-center justify-between gap-2";
    div.innerHTML = `
      <div>
        <div class="flex flex-wrap items-center gap-2">
          <span class="font-heading font-semibold text-sm text-slate-800">${idx + 1}. ${dv.label}</span>
          <span class="text-[11px] bg-emerald-100 text-emerald-800 px-2 py-0.5 rounded font-medium">5 ระดับ (1-5 คะแนน)</span>
        </div>
        <p class="text-xs text-slate-500 mt-1">ระดับคะแนน: 5=มากที่สุด, 4=มาก, 3=ปานกลาง, 2=น้อย, 1=น้อยที่สุด</p>
      </div>
      ${
        appState.dvs.length > 1
          ? `<button data-remove-dv="${dv.id}" class="text-xs text-rose-500 hover:text-rose-700 px-1.5 py-0.5" title="ลบตัวแปร">✕</button>`
          : ""
      }
    `;
    dvContainer.appendChild(div);
  });

  document.querySelectorAll("[data-remove-iv]").forEach((btn) => {
    btn.addEventListener("click", () => {
      const id = btn.getAttribute("data-remove-iv");
      appState.ivs = appState.ivs.filter((v) => v.id !== id);
      saveStateToLocal();
      renderAll();
    });
  });

  document.querySelectorAll("[data-remove-dv]").forEach((btn) => {
    btn.addEventListener("click", () => {
      const id = btn.getAttribute("data-remove-dv");
      appState.dvs = appState.dvs.filter((v) => v.id !== id);
      saveStateToLocal();
      renderAll();
    });
  });
}

function renderSurveyForm() {
  const ivFields = document.getElementById("formIvFields");
  ivFields.innerHTML = "";
  appState.ivs.forEach((iv) => {
    const wrapper = document.createElement("div");
    wrapper.innerHTML = `
      <label class="block text-xs font-medium text-slate-700 mb-1">${iv.label} <span class="text-sky-600">(ตัวแปรต้น)</span></label>
      <select name="iv_${iv.id}" required class="w-full border border-slate-300 rounded-lg px-3 py-2 text-xs sm:text-sm bg-white focus:outline-none focus:ring-2 focus:ring-brand-500">
        ${iv.options.map((opt) => `<option value="${opt}">${opt}</option>`).join("")}
      </select>
    `;
    ivFields.appendChild(wrapper);
  });

  // 5-Level Likert Radio Buttons for DVs
  const dvFields = document.getElementById("formDvFields");
  dvFields.innerHTML = "";
  appState.dvs.forEach((dv) => {
    const wrapper = document.createElement("div");
    wrapper.className = "bg-slate-50 p-3 rounded-xl border border-slate-200";
    wrapper.innerHTML = `
      <label class="block text-xs font-semibold text-slate-800 mb-2">${dv.label}</label>
      <div class="grid grid-cols-5 gap-1.5">
        ${[1, 2, 3, 4, 5]
          .map((lvl) => {
            const info = LIKERT_LEVELS.find((l) => l.val === lvl);
            return `
            <label class="likert-option cursor-pointer">
              <input type="radio" name="dv_${dv.id}" value="${lvl}" ${lvl === 4 ? "checked" : ""} class="sr-only" />
              <span class="flex flex-col items-center justify-center py-1.5 px-1 rounded-lg border border-slate-300 bg-white hover:bg-emerald-50 text-center transition">
                <strong class="text-sm font-heading">${lvl}</strong>
                <span class="text-[10px] leading-tight mt-0.5">${info.short}</span>
              </span>
            </label>
          `;
          })
          .join("")}
      </div>
    `;
    dvFields.appendChild(wrapper);
  });
}

function renderRawDataTable() {
  const thead = document.getElementById("rawDataThead");
  const tbody = document.getElementById("rawDataTbody");

  thead.innerHTML = `
    <tr>
      <th class="py-2.5 px-3 font-heading border-b border-slate-200 text-center w-14">คนที่</th>
      ${appState.ivs
        .map(
          (iv) =>
            `<th class="py-2.5 px-3 font-heading border-b border-slate-200 bg-sky-50/80 text-sky-950">
              <div class="text-[10px] text-sky-600 uppercase">ตัวแปรต้น (IV)</div>
              ${iv.label}
            </th>`
        )
        .join("")}
      ${appState.dvs
        .map(
          (dv) =>
            `<th class="py-2.5 px-3 font-heading border-b border-slate-200 bg-emerald-50/80 text-emerald-950">
              <div class="text-[10px] text-emerald-600 uppercase">ตัวแปรตาม (5 ระดับ)</div>
              ${dv.label}
            </th>`
        )
        .join("")}
      <th class="py-2.5 px-3 font-heading border-b border-slate-200 text-center w-12">ลบ</th>
    </tr>
  `;

  tbody.innerHTML = "";
  if (appState.data.length === 0) {
    tbody.innerHTML = `
      <tr>
        <td colspan="${appState.ivs.length + appState.dvs.length + 2}" class="py-8 text-center text-slate-400">
          ยังไม่มีข้อมูลในตาราง กรุณาเลือกตอบแบบสอบถามด้านซ้าย หรือกดปุ่ม "สุ่มเพิ่มข้อมูล"
        </td>
      </tr>
    `;
    return;
  }

  appState.data.forEach((row, rowIdx) => {
    const tr = document.createElement("tr");
    tr.className = "hover:bg-slate-50 transition";

    let cellsHtml = `<td class="py-2 px-3 text-center font-medium text-slate-500">${rowIdx + 1}</td>`;

    appState.ivs.forEach((iv) => {
      cellsHtml += `
        <td class="py-1.5 px-3">
          <select data-row="${rowIdx}" data-col="${iv.id}" data-type="iv" class="border border-transparent hover:border-slate-300 focus:border-brand-500 rounded px-1.5 py-1 bg-transparent text-xs sm:text-sm w-full">
            ${iv.options
              .map(
                (opt) =>
                  `<option value="${opt}" ${row[iv.id] === opt ? "selected" : ""}>${opt}</option>`
              )
              .join("")}
          </select>
        </td>
      `;
    });

    // 5-Level Select Dropdown in Data View
    appState.dvs.forEach((dv) => {
      const currentVal = clampLikert(row[dv.id]);
      cellsHtml += `
        <td class="py-1.5 px-3">
          <select data-row="${rowIdx}" data-col="${dv.id}" data-type="dv" class="border border-slate-200 hover:border-emerald-400 focus:border-emerald-600 rounded-lg px-2 py-1 bg-white text-xs sm:text-sm font-medium text-slate-800 w-full">
            ${LIKERT_LEVELS.map(
              (lvl) =>
                `<option value="${lvl.val}" ${currentVal === lvl.val ? "selected" : ""}>${lvl.label}</option>`
            ).join("")}
          </select>
        </td>
      `;
    });

    cellsHtml += `
      <td class="py-1.5 px-3 text-center">
        <button data-delete-row="${rowIdx}" class="text-slate-400 hover:text-rose-600 transition" title="ลบแถวนี้">🗑️</button>
      </td>
    `;

    tr.innerHTML = cellsHtml;
    tbody.appendChild(tr);
  });

  // Bind inline edits
  tbody.querySelectorAll("select[data-type='iv']").forEach((sel) => {
    sel.addEventListener("change", (e) => {
      const r = Number(sel.getAttribute("data-row"));
      const col = sel.getAttribute("data-col");
      appState.data[r][col] = e.target.value;
      saveStateToLocal();
      runAllAnalyses();
    });
  });

  tbody.querySelectorAll("select[data-type='dv']").forEach((sel) => {
    sel.addEventListener("change", (e) => {
      const r = Number(sel.getAttribute("data-row"));
      const col = sel.getAttribute("data-col");
      appState.data[r][col] = clampLikert(e.target.value);
      saveStateToLocal();
      runAllAnalyses();
    });
  });

  tbody.querySelectorAll("[data-delete-row]").forEach((btn) => {
    btn.addEventListener("click", () => {
      const r = Number(btn.getAttribute("data-delete-row"));
      appState.data.splice(r, 1);
      saveStateToLocal();
      renderProjectHeader();
      renderRawDataTable();
      runAllAnalyses();
    });
  });
}

// Simulate realistic random respondent with 5-level rating (1..5)
function simulateRespondent() {
  const newRow = {};
  appState.ivs.forEach((iv) => {
    const randIdx = Math.floor(Math.random() * iv.options.length);
    newRow[iv.id] = iv.options[randIdx];
  });

  appState.dvs.forEach((dv, idx) => {
    const existing = appState.data.map((r) => Number(r[dv.id])).filter((v) => !isNaN(v));
    const baseMean = existing.length > 1 ? StatMath.mean(existing) : 3.6;
    const firstIv = appState.ivs[0];
    const groupIdx = firstIv ? firstIv.options.indexOf(newRow[firstIv.id]) : 0;
    const u1 = Math.max(1e-7, Math.random());
    const u2 = Math.random();
    const z = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
    const rawScore = baseMean + z * 0.85 + groupIdx * (idx % 2 === 0 ? 0.45 : -0.35);
    newRow[dv.id] = clampLikert(rawScore);
  });

  appState.data.push(newRow);
}

// ============================================================================
// 6. POPULATE SELECTORS & RUN ALL 4 STATISTICAL ANALYSES
// ============================================================================

function populateAnalysisSelectors() {
  const freqSel = document.getElementById("freqVarSelect");
  const prevFreq = freqSel.value;
  freqSel.innerHTML = appState.ivs
    .map((iv) => `<option value="${iv.id}">${iv.label}</option>`)
    .join("");
  if (appState.ivs.some((v) => v.id === prevFreq)) freqSel.value = prevFreq;

  const oneDvSel = document.getElementById("oneSampleDvSelect");
  const prevOneDv = oneDvSel.value;
  oneDvSel.innerHTML = appState.dvs
    .map((dv) => `<option value="${dv.id}">${dv.label}</option>`)
    .join("");
  if (appState.dvs.some((v) => v.id === prevOneDv)) {
    oneDvSel.value = prevOneDv;
  } else if (appState.dvs[0]) {
    document.getElementById("oneSampleTestVal").value = appState.dvs[0].defaultTestVal ?? 3.5;
  }

  const indepIvSel = document.getElementById("indepIvSelect");
  const prevIndepIv = indepIvSel.value;
  indepIvSel.innerHTML = appState.ivs
    .map((iv) => `<option value="${iv.id}">${iv.label} (${iv.options.length} กลุ่ม)</option>`)
    .join("");
  if (appState.ivs.some((v) => v.id === prevIndepIv)) {
    indepIvSel.value = prevIndepIv;
  } else {
    const twoGroupIv = appState.ivs.find((v) => v.options.length === 2);
    if (twoGroupIv) indepIvSel.value = twoGroupIv.id;
  }
  updateIndepGroupSelectors();

  const indepDvSel = document.getElementById("indepDvSelect");
  const prevIndepDv = indepDvSel.value;
  indepDvSel.innerHTML = appState.dvs
    .map((dv) => `<option value="${dv.id}">${dv.label}</option>`)
    .join("");
  if (appState.dvs.some((v) => v.id === prevIndepDv)) indepDvSel.value = prevIndepDv;

  const anovaIvSel = document.getElementById("anovaIvSelect");
  const prevAnovaIv = anovaIvSel.value;
  anovaIvSel.innerHTML = appState.ivs
    .map((iv) => `<option value="${iv.id}">${iv.label} (${iv.options.length} กลุ่ม)</option>`)
    .join("");
  if (appState.ivs.some((v) => v.id === prevAnovaIv)) {
    anovaIvSel.value = prevAnovaIv;
  } else {
    const multiGroupIv = appState.ivs.find((v) => v.options.length >= 3);
    if (multiGroupIv) anovaIvSel.value = multiGroupIv.id;
  }

  const anovaDvSel = document.getElementById("anovaDvSelect");
  const prevAnovaDv = anovaDvSel.value;
  anovaDvSel.innerHTML = appState.dvs
    .map((dv) => `<option value="${dv.id}">${dv.label}</option>`)
    .join("");
  if (appState.dvs.some((v) => v.id === prevAnovaDv)) anovaDvSel.value = prevAnovaDv;
}

function updateIndepGroupSelectors() {
  const ivId = document.getElementById("indepIvSelect").value;
  const iv = appState.ivs.find((v) => v.id === ivId) || appState.ivs[0];
  if (!iv) return;

  const g1Sel = document.getElementById("indepGroup1Select");
  const g2Sel = document.getElementById("indepGroup2Select");
  g1Sel.innerHTML = iv.options.map((opt) => `<option value="${opt}">${opt}</option>`).join("");
  g2Sel.innerHTML = iv.options.map((opt) => `<option value="${opt}">${opt}</option>`).join("");

  if (iv.options.length >= 2) {
    g1Sel.value = iv.options[0];
    g2Sel.value = iv.options[1];
  }
}

function runAllAnalyses() {
  renderFrequencyAnalysis();
  renderOneSampleAnalysis();
  renderIndepTTestAnalysis();
  renderAnovaAnalysis();
}

// ----------------------------------------------------------------------------
// 6.1 Frequency & Percentage + 5-Level DV Distribution Table
// ----------------------------------------------------------------------------
function renderFrequencyAnalysis() {
  const ivId = document.getElementById("freqVarSelect").value;
  const iv = appState.ivs.find((v) => v.id === ivId) || appState.ivs[0];
  if (!iv) return;

  const N = appState.data.length;
  document.getElementById("freqTableTitle").textContent = `ตารางที่ 1 จำนวนและร้อยละของกลุ่มตัวอย่าง จำแนกตาม "${iv.label}" (N = ${N})`;

  let cumCount = 0;
  const rowsData = iv.options.map((opt) => {
    const count = appState.data.filter((r) => r[iv.id] === opt).length;
    cumCount += count;
    const pct = N > 0 ? (count / N) * 100 : 0;
    const cumPct = N > 0 ? (cumCount / N) * 100 : 0;
    return { option: opt, count, pct, cumPct };
  });

  const tbody = document.getElementById("freqTableBody");
  tbody.innerHTML = rowsData
    .map(
      (r) => `
      <tr class="hover:bg-slate-50">
        <td class="py-2.5 px-4 font-medium text-slate-800">${r.option}</td>
        <td class="py-2.5 px-4 text-right">${r.count}</td>
        <td class="py-2.5 px-4 text-right font-semibold text-indigo-700">${fmt(r.pct, 2)}</td>
        <td class="py-2.5 px-4 text-right text-slate-500">${fmt(r.cumPct, 2)}</td>
      </tr>
    `
    )
    .join("");

  document.getElementById("freqTableFoot").innerHTML = `
    <tr>
      <td class="py-2.5 px-4">รวมทั้งสิ้น (Total)</td>
      <td class="py-2.5 px-4 text-right">${N}</td>
      <td class="py-2.5 px-4 text-right">${N > 0 ? "100.00" : "0.00"}</td>
      <td class="py-2.5 px-4 text-right">-</td>
    </tr>
  `;

  // Bar Chart
  const chartContainer = document.getElementById("freqBarChart");
  chartContainer.innerHTML = rowsData
    .map(
      (r) => `
      <div>
        <div class="flex justify-between text-xs mb-1">
          <span class="font-medium text-slate-700">${r.option}</span>
          <span class="text-slate-600"><strong>${r.count}</strong> คน (${fmt(r.pct, 2)}%)</span>
        </div>
        <div class="w-full bg-slate-200 h-3 rounded-full overflow-hidden">
          <div class="bg-brand-600 h-full rounded-full transition-all duration-300" style="width: ${Math.max(2, r.pct)}%"></div>
        </div>
      </div>
    `
    )
    .join("");

  // Table 2: 5-Level DV Distribution Table (Count & % for levels 5, 4, 3, 2, 1 + Mean, SD, Meaning)
  const descDvBody = document.getElementById("descDvTableBody");
  descDvBody.innerHTML = appState.dvs
    .map((dv) => {
      const vals = appState.data.map((r) => clampLikert(r[dv.id]));
      const n = vals.length;
      const mean = StatMath.mean(vals);
      const sd = StatMath.sd(vals);
      const meaning = interpretMean5Level(mean);

      const levelCells = [5, 4, 3, 2, 1]
        .map((lvl) => {
          const c = vals.filter((v) => v === lvl).length;
          const pct = n > 0 ? (c / n) * 100 : 0;
          return `<td class="py-2 px-2 text-center border-r border-slate-200 text-xs">
            <span class="font-semibold text-slate-800">${c}</span>
            <span class="text-slate-500">(${fmt(pct, 1)}%)</span>
          </td>`;
        })
        .join("");

      return `
        <tr class="hover:bg-slate-50">
          <td class="py-2.5 px-3 font-medium text-slate-800 border-r border-slate-300">${dv.label}</td>
          ${levelCells}
          <td class="py-2.5 px-3 text-right font-bold text-emerald-700">${fmt(mean, 2)}</td>
          <td class="py-2.5 px-3 text-right">${fmt(sd, 2)}</td>
          <td class="py-2.5 px-3 text-center">${meaning.badge}</td>
        </tr>
      `;
    })
    .join("");

  // Interpretation Text
  const sorted = [...rowsData].sort((a, b) => b.count - a.count);
  if (N === 0) {
    document.getElementById("freqInterpretationText").textContent = "ยังไม่มีข้อมูลสำหรับการแปลผล";
  } else {
    const parts = sorted.map(
      (r, idx) =>
        `${idx === 0 ? "ส่วนใหญ่เป็นกลุ่ม" : "รองลงมาคือกลุ่ม"} <strong>${r.option}</strong> จำนวน ${r.count} คน (ร้อยละ ${fmt(r.pct, 2)})`
    );
    document.getElementById("freqInterpretationText").innerHTML =
      `จากตารางที่ 1 พบว่า ผู้ตอบแบบสอบถามทั้งหมดจำนวน ${N} คน เมื่อจำแนกตาม<strong>${iv.label}</strong> ${parts.join(" และ ")} ตามลำดับ`;
  }
}

// ----------------------------------------------------------------------------
// 6.2 One-Sample t-test (5-Level DV vs Criterion Test Value)
// ----------------------------------------------------------------------------
function renderOneSampleAnalysis() {
  const dvId = document.getElementById("oneSampleDvSelect").value;
  const dv = appState.dvs.find((v) => v.id === dvId) || appState.dvs[0];
  if (!dv) return;

  const mu0 = parseFloat(document.getElementById("oneSampleTestVal").value) || 3.5;
  const alpha = parseFloat(document.getElementById("oneSampleAlpha").value) || 0.05;

  document.getElementById("oneSampleTestValBadge").textContent = `Test Value (เกณฑ์เปรียบเทียบ μ₀) = ${fmt(mu0, 2)} คะแนน`;
  document.getElementById("oneSampleH0").innerHTML = `<strong>H₀ : μ = ${fmt(mu0, 2)}</strong> (ค่าเฉลี่ยของ${dv.label} <strong>ไม่แตกต่าง</strong>จากเกณฑ์ ${fmt(mu0, 2)} คะแนน)`;
  document.getElementById("oneSampleH1").innerHTML = `<strong>H₁ : μ ≠ ${fmt(mu0, 2)}</strong> (ค่าเฉลี่ยของ${dv.label} <strong>แตกต่าง</strong>จากเกณฑ์ ${fmt(mu0, 2)} คะแนน)`;

  const vals = appState.data.map((r) => clampLikert(r[dv.id]));
  const res = StatMath.oneSampleTTest(vals, mu0, alpha);
  const meaning = interpretMean5Level(res.mean);

  document.getElementById("oneSampleStatBody").innerHTML = `
    <tr>
      <td class="py-2.5 px-4 font-medium text-slate-800">${dv.label} ${meaning.badge}</td>
      <td class="py-2.5 px-4 text-right">${res.n}</td>
      <td class="py-2.5 px-4 text-right font-semibold text-indigo-700">${fmt(res.mean, 2)}</td>
      <td class="py-2.5 px-4 text-right">${fmt(res.sd, 2)}</td>
      <td class="py-2.5 px-4 text-right text-slate-500">${fmt(res.se, 3)}</td>
    </tr>
  `;

  const isSig = res.p < alpha;
  document.getElementById("oneSampleTestBody").innerHTML = `
    <tr class="${isSig ? "bg-amber-50/70" : ""}">
      <td class="py-2.5 px-4 font-medium text-slate-800">${dv.label}</td>
      <td class="py-2.5 px-4 text-right font-semibold">${fmt(res.t, 3)}</td>
      <td class="py-2.5 px-4 text-right">${res.df}</td>
      <td class="py-2.5 px-4 text-right font-bold ${isSig ? "text-emerald-700" : "text-slate-700"}">${fmtP(res.p)}</td>
      <td class="py-2.5 px-4 text-right">${fmt(res.diff, 2)}</td>
      <td class="py-2.5 px-4 text-right text-slate-600">[${fmt(res.ciLow, 2)}, ${fmt(res.ciHigh, 2)}]</td>
    </tr>
  `;

  const badge = document.getElementById("oneSampleDecisionBadge");
  if (isSig) {
    badge.className = "text-xs font-bold px-2.5 py-0.5 rounded-full bg-emerald-100 text-emerald-800 border border-emerald-300";
    badge.textContent = `Sig. < ${alpha} → ปฏิเสธ H₀ (แตกต่างอย่างมีนัยสำคัญ)`;
  } else {
    badge.className = "text-xs font-bold px-2.5 py-0.5 rounded-full bg-slate-200 text-slate-700";
    badge.textContent = `Sig. ≥ ${alpha} → ยอมรับ H₀ (ไม่แตกต่างจากเกณฑ์)`;
  }

  const dirText = res.diff >= 0 ? "สูงกว่า" : "ต่ำกว่า";
  document.getElementById("oneSampleSummaryText").innerHTML = isSig
    ? `ผลการวิเคราะห์ด้วยสถิติ One-Sample t-test พบว่า ผู้ตอบแบบสอบถาม (${res.n} คน) มี<strong>${dv.label}</strong>อยู่ในระดับ<strong>${meaning.text}</strong> (<span class="xbar">x</span> = ${fmt(res.mean, 2)}, S.D. = ${fmt(res.sd, 2)}) เมื่อเปรียบเทียบกับเกณฑ์ ${fmt(mu0, 2)} คะแนน พบว่า มีค่าเฉลี่ย<strong>${dirText}เกณฑ์ที่กำหนดอย่างมีนัยสำคัญทางสถิติที่ระดับ ${alpha}</strong> (<i>t</i> = ${fmt(res.t, 3)}, <i>df</i> = ${res.df}, <i>p</i> = ${fmtPRaw(res.p)})`
    : `ผลการวิเคราะห์ด้วยสถิติ One-Sample t-test พบว่า ผู้ตอบแบบสอบถาม (${res.n} คน) มี<strong>${dv.label}</strong>อยู่ในระดับ<strong>${meaning.text}</strong> (<span class="xbar">x</span> = ${fmt(res.mean, 2)}, S.D. = ${fmt(res.sd, 2)}) เมื่อเปรียบเทียบกับเกณฑ์ ${fmt(mu0, 2)} คะแนน พบว่า <strong>ไม่แตกต่างจากเกณฑ์ที่กำหนด</strong>ที่ระดับนัยสำคัญ ${alpha} (<i>t</i> = ${fmt(res.t, 3)}, <i>df</i> = ${res.df}, <i>p</i> = ${fmtPRaw(res.p)})`;
}

// ----------------------------------------------------------------------------
// 6.3 Independent Samples t-test (Full SPSS 2-row table + APA table)
// ----------------------------------------------------------------------------
function renderIndepTTestAnalysis() {
  const ivId = document.getElementById("indepIvSelect").value;
  const dvId = document.getElementById("indepDvSelect").value;
  const iv = appState.ivs.find((v) => v.id === ivId) || appState.ivs[0];
  const dv = appState.dvs.find((v) => v.id === dvId) || appState.dvs[0];
  if (!iv || !dv) return;

  const g1Name = document.getElementById("indepGroup1Select").value;
  const g2Name = document.getElementById("indepGroup2Select").value;

  document.getElementById("indepH0").innerHTML = `<strong>H₀ : μ₁ = μ₂</strong> (กลุ่ม<strong>${g1Name}</strong> และ <strong>${g2Name}</strong> มีค่าเฉลี่ย${dv.label}<strong>ไม่แตกต่างกัน</strong>)`;
  document.getElementById("indepH1").innerHTML = `<strong>H₁ : μ₁ ≠ μ₂</strong> (กลุ่ม<strong>${g1Name}</strong> และ <strong>${g2Name}</strong> มีค่าเฉลี่ย${dv.label}<strong>แตกต่างกัน</strong>)`;

  const g1Vals = appState.data
    .filter((r) => r[iv.id] === g1Name)
    .map((r) => clampLikert(r[dv.id]));
  const g2Vals = appState.data
    .filter((r) => r[iv.id] === g2Name)
    .map((r) => clampLikert(r[dv.id]));

  const res = StatMath.independentTTest(g1Vals, g2Vals);
  const m1Info = interpretMean5Level(res.g1.mean);
  const m2Info = interpretMean5Level(res.g2.mean);

  // Table 1: Group Statistics
  document.getElementById("indepGroupStatBody").innerHTML = `
    <tr>
      <td rowspan="2" class="py-2.5 px-4 font-medium text-slate-800 border-r border-slate-200">${dv.label}</td>
      <td class="py-2.5 px-4 font-semibold text-sky-900">${g1Name}</td>
      <td class="py-2.5 px-4 text-right">${res.g1.n}</td>
      <td class="py-2.5 px-4 text-right font-semibold text-indigo-700">${fmt(res.g1.mean, 2)}</td>
      <td class="py-2.5 px-4 text-right">${fmt(res.g1.sd, 2)}</td>
      <td class="py-2.5 px-4 text-center">${m1Info.badge}</td>
    </tr>
    <tr>
      <td class="py-2.5 px-4 font-semibold text-sky-900">${g2Name}</td>
      <td class="py-2.5 px-4 text-right">${res.g2.n}</td>
      <td class="py-2.5 px-4 text-right font-semibold text-indigo-700">${fmt(res.g2.mean, 2)}</td>
      <td class="py-2.5 px-4 text-right">${fmt(res.g2.sd, 2)}</td>
      <td class="py-2.5 px-4 text-center">${m2Info.badge}</td>
    </tr>
  `;

  // Table 2: SPSS 2-Row Independent Samples Test
  const eqClass = res.equalAssumed ? "selected-levene-row" : "text-slate-500";
  const neqClass = !res.equalAssumed ? "selected-levene-row" : "text-slate-500";

  document.getElementById("indepSpssTestBody").innerHTML = `
    <tr class="${eqClass}">
      <td class="py-2.5 px-3 border-r border-slate-300">
        Equal variances assumed (ความแปรปรวนเท่ากัน) ${res.equalAssumed ? "👈 <span class='text-xs text-amber-800'>(อ่านบรรทัดนี้)</span>" : ""}
      </td>
      <td class="py-2.5 px-3 text-right bg-sky-50/40 font-semibold">${fmt(res.leveneF, 3)}</td>
      <td class="py-2.5 px-3 text-right border-r border-slate-300 bg-sky-50/40 font-semibold">${fmtPRaw(res.leveneP)}</td>
      <td class="py-2.5 px-3 text-right">${fmt(res.equalRow.t, 3)}</td>
      <td class="py-2.5 px-3 text-right">${res.equalRow.df}</td>
      <td class="py-2.5 px-3 text-right font-bold">${fmtP(res.equalRow.p)}</td>
      <td class="py-2.5 px-3 text-right">${fmt(res.meanDiff, 2)}</td>
      <td class="py-2.5 px-3 text-right">${fmt(res.equalRow.seDiff, 3)}</td>
    </tr>
    <tr class="${neqClass}">
      <td class="py-2.5 px-3 border-r border-slate-300">
        Equal variances not assumed (ความแปรปรวนไม่เท่ากัน) ${!res.equalAssumed ? "👈 <span class='text-xs text-amber-800'>(อ่านบรรทัดนี้)</span>" : ""}
      </td>
      <td class="py-2.5 px-3 text-right bg-sky-50/20 text-slate-300">-</td>
      <td class="py-2.5 px-3 text-right border-r border-slate-300 bg-sky-50/20 text-slate-300">-</td>
      <td class="py-2.5 px-3 text-right">${fmt(res.unequalRow.t, 3)}</td>
      <td class="py-2.5 px-3 text-right">${fmt(res.unequalRow.df, 2)}</td>
      <td class="py-2.5 px-3 text-right font-bold">${fmtP(res.unequalRow.p)}</td>
      <td class="py-2.5 px-3 text-right">${fmt(res.meanDiff, 2)}</td>
      <td class="py-2.5 px-3 text-right">${fmt(res.unequalRow.seDiff, 3)}</td>
    </tr>
  `;

  const chosen = res.equalAssumed ? res.equalRow : res.unequalRow;
  const chosenDfStr = res.equalAssumed ? chosen.df : fmt(chosen.df, 2);

  // Table 3: APA Summary Table
  document.getElementById("indepApaTableBody").innerHTML = `
    <tr>
      <td class="py-2.5 px-4 font-medium">${g1Name}</td>
      <td class="py-2.5 px-4 text-right">${res.g1.n}</td>
      <td class="py-2.5 px-4 text-right font-semibold">${fmt(res.g1.mean, 2)}</td>
      <td class="py-2.5 px-4 text-right">${fmt(res.g1.sd, 2)}</td>
      <td rowspan="2" class="py-2.5 px-4 text-right font-semibold align-middle border-l border-slate-200">${fmt(chosen.t, 3)}</td>
      <td rowspan="2" class="py-2.5 px-4 text-right align-middle">${chosenDfStr}</td>
      <td rowspan="2" class="py-2.5 px-4 text-right font-bold align-middle ${chosen.p < 0.05 ? "text-emerald-700" : ""}">${fmtP(chosen.p)}</td>
    </tr>
    <tr>
      <td class="py-2.5 px-4 font-medium">${g2Name}</td>
      <td class="py-2.5 px-4 text-right">${res.g2.n}</td>
      <td class="py-2.5 px-4 text-right font-semibold">${fmt(res.g2.mean, 2)}</td>
      <td class="py-2.5 px-4 text-right">${fmt(res.g2.sd, 2)}</td>
    </tr>
  `;

  const isSig = chosen.p < 0.05;
  const badge = document.getElementById("indepDecisionBadge");
  badge.className = isSig
    ? "text-xs font-bold px-2.5 py-0.5 rounded-full bg-emerald-100 text-emerald-800 border border-emerald-300"
    : "text-xs font-bold px-2.5 py-0.5 rounded-full bg-slate-200 text-slate-700";
  badge.textContent = isSig ? "Sig. < 0.05 → แตกต่างกันอย่างมีนัยสำคัญ" : "Sig. ≥ 0.05 → ไม่แตกต่างกัน";

  document.getElementById("indepStepWalkthrough").innerHTML = `
    <div><strong>ขั้นที่ 1 (ตรวจ Levene's Test):</strong> F = ${fmt(res.leveneF, 3)}, Sig. = ${fmtPRaw(res.leveneP)} ${
    res.equalAssumed
      ? "(มากกว่า 0.05) ดังนั้นความแปรปรวนเท่ากัน → <strong>เลือกอ่านบรรทัดที่ 1 (Equal variances assumed)</strong>"
      : "(น้อยกว่าหรือเท่ากับ 0.05) ดังนั้นความแปรปรวนไม่เท่ากัน → <strong>เลือกอ่านบรรทัดที่ 2 (Equal variances not assumed)</strong>"
  }</div>
    <div><strong>ขั้นที่ 2 (อ่านผล t-test บรรทัดที่เลือก):</strong> <i>t</i> = ${fmt(chosen.t, 3)}, <i>df</i> = ${chosenDfStr}, Sig. (2-tailed) = ${fmtPRaw(chosen.p)}</div>
  `;

  const higherGroup = res.g1.mean >= res.g2.mean ? g1Name : g2Name;
  const lowerGroup = res.g1.mean >= res.g2.mean ? g2Name : g1Name;

  document.getElementById("indepSummaryText").innerHTML = isSig
    ? `ผลการเปรียบเทียบ<strong>${dv.label}</strong> จำแนกตาม<strong>${iv.label}</strong> ระหว่างกลุ่ม<strong>${g1Name}</strong> (<span class="xbar">x</span> = ${fmt(res.g1.mean, 2)}, S.D. = ${fmt(res.g1.sd, 2)}) และกลุ่ม<strong>${g2Name}</strong> (<span class="xbar">x</span> = ${fmt(res.g2.mean, 2)}, S.D. = ${fmt(res.g2.sd, 2)}) พบว่า <strong>แตกต่างกันอย่างมีนัยสำคัญทางสถิติที่ระดับ .05</strong> (<i>t</i> = ${fmt(chosen.t, 3)}, <i>df</i> = ${chosenDfStr}, <i>p</i> = ${fmtPRaw(chosen.p)}) โดยกลุ่ม<strong>${higherGroup}</strong>มีค่าเฉลี่ยสูงกว่ากลุ่ม<strong>${lowerGroup}</strong>`
    : `ผลการเปรียบเทียบ<strong>${dv.label}</strong> จำแนกตาม<strong>${iv.label}</strong> ระหว่างกลุ่ม<strong>${g1Name}</strong> (<span class="xbar">x</span> = ${fmt(res.g1.mean, 2)}, S.D. = ${fmt(res.g1.sd, 2)}) และกลุ่ม<strong>${g2Name}</strong> (<span class="xbar">x</span> = ${fmt(res.g2.mean, 2)}, S.D. = ${fmt(res.g2.sd, 2)}) พบว่า <strong>ไม่แตกต่างกัน</strong>ทางสถิติ (<i>t</i> = ${fmt(chosen.t, 3)}, <i>df</i> = ${chosenDfStr}, <i>p</i> = ${fmtPRaw(chosen.p)})`;
}

// ----------------------------------------------------------------------------
// 6.4 One-Way ANOVA (Descriptives + ANOVA Table + Post-hoc LSD)
// ----------------------------------------------------------------------------
function renderAnovaAnalysis() {
  const ivId = document.getElementById("anovaIvSelect").value;
  const dvId = document.getElementById("anovaDvSelect").value;
  const iv = appState.ivs.find((v) => v.id === ivId) || appState.ivs[0];
  const dv = appState.dvs.find((v) => v.id === dvId) || appState.dvs[0];
  if (!iv || !dv) return;

  document.getElementById("anovaH0").innerHTML = `<strong>H₀ : μ₁ = μ₂ = … = μₖ</strong> (ค่าเฉลี่ยของ${dv.label} จำแนกตาม${iv.label}ทุกกลุ่ม <strong>ไม่แตกต่างกัน</strong>)`;
  document.getElementById("anovaH1").innerHTML = `<strong>H₁ :</strong> มีอย่างน้อย 1 คู่กลุ่มที่มีค่าเฉลี่ยของ${dv.label} <strong>แตกต่างกัน</strong>`;

  const groupsObj = iv.options.map((opt) => ({
    name: opt,
    values: appState.data
      .filter((r) => r[iv.id] === opt)
      .map((r) => clampLikert(r[dv.id]))
  }));

  const res = StatMath.oneWayAnova(groupsObj);
  const grandMeaning = interpretMean5Level(res.grandMean);

  // Table 1: Descriptives
  document.getElementById("anovaDescBody").innerHTML =
    res.descriptives
      .map((d) => {
        const mInfo = interpretMean5Level(d.mean);
        return `
        <tr class="hover:bg-slate-50">
          <td class="py-2.5 px-4 font-semibold text-sky-900">${d.name}</td>
          <td class="py-2.5 px-4 text-right">${d.n}</td>
          <td class="py-2.5 px-4 text-right font-semibold text-indigo-700">${fmt(d.mean, 2)}</td>
          <td class="py-2.5 px-4 text-right">${fmt(d.sd, 2)}</td>
          <td class="py-2.5 px-4 text-center">${mInfo.badge}</td>
          <td class="py-2.5 px-4 text-right">${d.min}</td>
          <td class="py-2.5 px-4 text-right">${d.max}</td>
        </tr>
      `;
      })
      .join("") +
    `
      <tr class="bg-slate-50 font-bold border-t-2 border-slate-800">
        <td class="py-2.5 px-4">รวมทั้งหมด (Total)</td>
        <td class="py-2.5 px-4 text-right">${res.N}</td>
        <td class="py-2.5 px-4 text-right text-indigo-800">${fmt(res.grandMean, 2)}</td>
        <td class="py-2.5 px-4 text-right">${fmt(res.grandSd, 2)}</td>
        <td class="py-2.5 px-4 text-center">${grandMeaning.badge}</td>
        <td class="py-2.5 px-4 text-right">-</td>
        <td class="py-2.5 px-4 text-right">-</td>
      </tr>
    `;

  const isSig = res.pVal < 0.05;

  // Table 2: ANOVA Table
  document.getElementById("anovaTableBody").innerHTML = `
    <tr class="${isSig ? "bg-amber-50/70" : ""}">
      <td class="py-2.5 px-4 font-semibold text-slate-800">ระหว่างกลุ่ม (Between Groups)</td>
      <td class="py-2.5 px-4 text-right">${fmt(res.ssBetween, 3)}</td>
      <td class="py-2.5 px-4 text-right">${res.dfBetween}</td>
      <td class="py-2.5 px-4 text-right">${fmt(res.msBetween, 3)}</td>
      <td class="py-2.5 px-4 text-right font-bold text-indigo-700">${fmt(res.fVal, 3)}</td>
      <td class="py-2.5 px-4 text-right font-bold ${isSig ? "text-emerald-700" : "text-slate-700"}">${fmtP(res.pVal)}</td>
    </tr>
    <tr>
      <td class="py-2.5 px-4 font-semibold text-slate-800">ภายในกลุ่ม (Within Groups / Error)</td>
      <td class="py-2.5 px-4 text-right">${fmt(res.ssWithin, 3)}</td>
      <td class="py-2.5 px-4 text-right">${res.dfWithin}</td>
      <td class="py-2.5 px-4 text-right">${fmt(res.msWithin, 3)}</td>
      <td class="py-2.5 px-4 text-right text-slate-400">-</td>
      <td class="py-2.5 px-4 text-right text-slate-400">-</td>
    </tr>
    <tr class="bg-slate-50 font-bold border-t-2 border-slate-800">
      <td class="py-2.5 px-4">รวม (Total)</td>
      <td class="py-2.5 px-4 text-right">${fmt(res.ssTotal, 3)}</td>
      <td class="py-2.5 px-4 text-right">${res.dfTotal}</td>
      <td class="py-2.5 px-4 text-right text-slate-400">-</td>
      <td class="py-2.5 px-4 text-right text-slate-400">-</td>
      <td class="py-2.5 px-4 text-right text-slate-400">-</td>
    </tr>
  `;

  // Table 3: Post-hoc Pairwise Comparison
  document.getElementById("anovaPosthocBody").innerHTML = res.posthoc
    .map((pair) => {
      const pairSig = pair.p < 0.05;
      return `
        <tr class="${pairSig && isSig ? "bg-emerald-50/60" : ""}">
          <td class="py-2.5 px-4 font-medium text-slate-800">${pair.groupI} vs ${pair.groupJ}</td>
          <td class="py-2.5 px-4 text-right font-semibold">${fmt(pair.diff, 2)}${pairSig ? "*" : ""}</td>
          <td class="py-2.5 px-4 text-right text-slate-500">${fmt(pair.se, 3)}</td>
          <td class="py-2.5 px-4 text-right font-bold ${pairSig ? "text-emerald-700" : "text-slate-600"}">${fmtP(pair.p)}</td>
          <td class="py-2.5 px-4 text-xs">
            ${
              pairSig
                ? `<span class="text-emerald-800 font-semibold">แตกต่างกันอย่างมีนัยสำคัญ</span> (${
                    pair.diff > 0 ? `${pair.groupI} > ${pair.groupJ}` : `${pair.groupJ} > ${pair.groupI}`
                  })`
                : `<span class="text-slate-500">ไม่แตกต่างกัน</span>`
            }
          </td>
        </tr>
      `;
    })
    .join("");

  const badge = document.getElementById("anovaDecisionBadge");
  badge.className = isSig
    ? "text-xs font-bold px-2.5 py-0.5 rounded-full bg-emerald-100 text-emerald-800 border border-emerald-300"
    : "text-xs font-bold px-2.5 py-0.5 rounded-full bg-slate-200 text-slate-700";
  badge.textContent = isSig
    ? "Sig. < 0.05 → ปฏิเสธ H₀ (มีอย่างน้อย 1 คู่ต่างกัน)"
    : "Sig. ≥ 0.05 → ยอมรับ H₀ (ไม่แตกต่างกัน)";

  const sigPairs = res.posthoc.filter((p) => p.p < 0.05);
  const sigPairsText =
    sigPairs.length > 0
      ? `เมื่อทดสอบรายคู่ด้วยวิธี LSD พบคู่ที่มีค่าเฉลี่ยแตกต่างกันอย่างมีนัยสำคัญทางสถิติที่ระดับ .05 จำนวน ${sigPairs.length} คู่ ได้แก่ ${sigPairs
          .map((p) => `คู่<strong>${p.groupI} กับ ${p.groupJ}</strong>`)
          .join(", ")}`
      : "";

  document.getElementById("anovaSummaryText").innerHTML = isSig
    ? `ผลการวิเคราะห์ความแปรปรวนทางเดียว (One-Way ANOVA) ของ<strong>${dv.label}</strong> จำแนกตาม<strong>${iv.label}</strong> (${res.k} กลุ่ม) พบว่า <strong>มีค่าเฉลี่ยแตกต่างกันอย่างมีนัยสำคัญทางสถิติที่ระดับ .05</strong> (<i>F</i> = ${fmt(res.fVal, 3)}, <i>df</i> = (${res.dfBetween}, ${res.dfWithin}), <i>p</i> = ${fmtPRaw(res.pVal)}) ${sigPairsText}`
    : `ผลการวิเคราะห์ความแปรปรวนทางเดียว (One-Way ANOVA) ของ<strong>${dv.label}</strong> จำแนกตาม<strong>${iv.label}</strong> (${res.k} กลุ่ม) พบว่า <strong>ไม่แตกต่างกัน</strong>ทางสถิติ (<i>F</i> = ${fmt(res.fVal, 3)}, <i>df</i> = (${res.dfBetween}, ${res.dfWithin}), <i>p</i> = ${fmtPRaw(res.pVal)}) จึงไม่จำเป็นต้องทดสอบรายคู่`;
}

// ============================================================================
// 7. INTERACTIVE PRACTICE QUIZ FOR STUDENTS (TAB 4)
// ============================================================================

function generateQuiz() {
  const topic = document.getElementById("quizTopicSelect").value;
  const feedbackBox = document.getElementById("quizFeedbackBox");
  feedbackBox.classList.add("hidden");
  feedbackBox.innerHTML = "";

  const tableContainer = document.getElementById("quizTableContainer");
  const qContainer = document.getElementById("quizQuestionsContainer");

  if (topic === "indep") {
    const iv = appState.ivs.find((v) => v.options.length === 2) || appState.ivs[0];
    const dv = appState.dvs[Math.floor(Math.random() * appState.dvs.length)];
    const g1 = iv.options[0];
    const g2 = iv.options[1] || iv.options[0];
    const g1Vals = appState.data.filter((r) => r[iv.id] === g1).map((r) => clampLikert(r[dv.id]));
    const g2Vals = appState.data.filter((r) => r[iv.id] === g2).map((r) => clampLikert(r[dv.id]));
    const res = StatMath.independentTTest(g1Vals, g2Vals);
    const chosen = res.equalAssumed ? res.equalRow : res.unequalRow;
    const isSig = chosen.p < 0.05;

    document.getElementById("quizScenarioTag").textContent = "โจทย์ฝึกอ่านตาราง: Independent Samples t-test";
    document.getElementById("quizScenarioTitle").textContent = `เปรียบเทียบ "${dv.label} (5 ระดับ)" ระหว่างกลุ่ม "${g1}" และ "${g2}"`;

    tableContainer.innerHTML = `
      <div class="text-xs font-semibold text-slate-700 mb-1">ตารางที่ 1 Group Statistics</div>
      <table class="w-full text-left border-collapse text-xs bg-white border border-slate-300 mb-3">
        <thead>
          <tr class="bg-slate-100 border-b border-slate-400">
            <th class="p-2">กลุ่ม</th><th class="p-2 text-right">n</th><th class="p-2 text-right"><span class="xbar">x</span> (Mean)</th><th class="p-2 text-right">S.D.</th>
          </tr>
        </thead>
        <tbody>
          <tr class="border-b border-slate-200"><td class="p-2 font-medium">${g1}</td><td class="p-2 text-right">${res.g1.n}</td><td class="p-2 text-right">${fmt(res.g1.mean, 2)}</td><td class="p-2 text-right">${fmt(res.g1.sd, 2)}</td></tr>
          <tr><td class="p-2 font-medium">${g2}</td><td class="p-2 text-right">${res.g2.n}</td><td class="p-2 text-right">${fmt(res.g2.mean, 2)}</td><td class="p-2 text-right">${fmt(res.g2.sd, 2)}</td></tr>
        </tbody>
      </table>

      <div class="text-xs font-semibold text-slate-700 mb-1">ตารางที่ 2 Independent Samples Test</div>
      <table class="w-full text-left border-collapse text-xs bg-white border border-slate-300">
        <thead>
          <tr class="bg-slate-100 border-b border-slate-400">
            <th class="p-2">Assumptions</th>
            <th class="p-2 text-right bg-sky-50">Levene's F</th>
            <th class="p-2 text-right bg-sky-50">Sig. (Levene)</th>
            <th class="p-2 text-right">t</th>
            <th class="p-2 text-right">df</th>
            <th class="p-2 text-right">Sig. (2-tailed)</th>
          </tr>
        </thead>
        <tbody>
          <tr class="border-b border-slate-200">
            <td class="p-2">Equal variances assumed</td>
            <td class="p-2 text-right bg-sky-50/40">${fmt(res.leveneF, 3)}</td>
            <td class="p-2 text-right bg-sky-50/40">${fmtPRaw(res.leveneP)}</td>
            <td class="p-2 text-right">${fmt(res.equalRow.t, 3)}</td>
            <td class="p-2 text-right">${res.equalRow.df}</td>
            <td class="p-2 text-right">${fmtPRaw(res.equalRow.p)}</td>
          </tr>
          <tr>
            <td class="p-2">Equal variances not assumed</td>
            <td class="p-2 text-right">-</td>
            <td class="p-2 text-right">-</td>
            <td class="p-2 text-right">${fmt(res.unequalRow.t, 3)}</td>
            <td class="p-2 text-right">${fmt(res.unequalRow.df, 2)}</td>
            <td class="p-2 text-right">${fmtPRaw(res.unequalRow.p)}</td>
          </tr>
        </tbody>
      </table>
    `;

    appState.currentQuiz = [
      {
        q: "1. จากค่า Sig. ของ Levene's Test นักศึกษาต้องเลือกอ่านค่า t จากบรรทัดใด?",
        options: [
          "บรรทัดบน (Equal variances assumed) เพราะ Sig. (Levene) > 0.05",
          "บรรทัดล่าง (Equal variances not assumed) เพราะ Sig. (Levene) ≤ 0.05"
        ],
        correct: res.equalAssumed ? 0 : 1,
        explanation: `ค่า Sig. ของ Levene's Test เท่ากับ ${fmtPRaw(res.leveneP)} ซึ่ง${res.equalAssumed ? "มากกว่า 0.05 แปลว่าความแปรปรวนเท่ากัน จึงต้องอ่านบรรทัดบน" : "น้อยกว่าหรือเท่ากับ 0.05 แปลว่าความแปรปรวนไม่เท่ากัน จึงต้องอ่านบรรทัดล่าง"}`
      },
      {
        q: `2. ค่าสถิติ t และองศาอิสระ (df) ที่ถูกต้องสำหรับการรายงานผลคือข้อใด?`,
        options: [
          `t = ${fmt(res.equalRow.t, 3)}, df = ${res.equalRow.df}`,
          `t = ${fmt(res.unequalRow.t, 3)}, df = ${fmt(res.unequalRow.df, 2)}`,
          `t = ${fmt(res.leveneF, 3)}, df = ${res.g1.n + res.g2.n}`
        ],
        correct: res.equalAssumed ? 0 : 1,
        explanation: `เมื่อเลือกบรรทัดที่ถูกต้อง จะได้ค่า t = ${fmt(chosen.t, 3)} และ df = ${res.equalAssumed ? chosen.df : fmt(chosen.df, 2)}`
      },
      {
        q: "3. ข้อใดสรุปผลการเปรียบเทียบที่ระดับนัยสำคัญ .05 ได้ถูกต้อง?",
        options: [
          `แตกต่างกันอย่างมีนัยสำคัญทางสถิติที่ระดับ .05 (เพราะ Sig. 2-tailed < 0.05)`,
          `ไม่แตกต่างกันทางสถิติ (เพราะ Sig. 2-tailed ≥ 0.05)`
        ],
        correct: isSig ? 0 : 1,
        explanation: `ค่า Sig. (2-tailed) เท่ากับ ${fmtPRaw(chosen.p)} ซึ่ง${isSig ? "น้อยกว่า 0.05 จึงปฏิเสธ H₀ (แตกต่างกันอย่างมีนัยสำคัญ)" : "มากกว่าหรือเท่ากับ 0.05 จึงยอมรับ H₀ (ไม่แตกต่างกัน)"}`
      }
    ];
  } else if (topic === "anova") {
    const iv = appState.ivs.find((v) => v.options.length >= 3) || appState.ivs[0];
    const dv = appState.dvs[Math.floor(Math.random() * appState.dvs.length)];
    const groupsObj = iv.options.map((opt) => ({
      name: opt,
      values: appState.data.filter((r) => r[iv.id] === opt).map((r) => clampLikert(r[dv.id]))
    }));
    const res = StatMath.oneWayAnova(groupsObj);
    const isSig = res.pVal < 0.05;

    document.getElementById("quizScenarioTag").textContent = "โจทย์ฝึกอ่านตาราง: One-Way ANOVA";
    document.getElementById("quizScenarioTitle").textContent = `เปรียบเทียบ "${dv.label} (5 ระดับ)" จำแนกตาม "${iv.label}" (${res.k} กลุ่ม)`;

    tableContainer.innerHTML = `
      <div class="text-xs font-semibold text-slate-700 mb-1">ตารางวิเคราะห์ความแปรปรวน (ANOVA Table)</div>
      <table class="w-full text-left border-collapse text-xs bg-white border border-slate-300">
        <thead>
          <tr class="bg-slate-100 border-b border-slate-400">
            <th class="p-2">Source</th><th class="p-2 text-right">SS</th><th class="p-2 text-right">df</th><th class="p-2 text-right">MS</th><th class="p-2 text-right">F</th><th class="p-2 text-right">Sig.</th>
          </tr>
        </thead>
        <tbody>
          <tr class="border-b border-slate-200">
            <td class="p-2 font-medium">Between Groups</td>
            <td class="p-2 text-right">${fmt(res.ssBetween, 3)}</td>
            <td class="p-2 text-right">${res.dfBetween}</td>
            <td class="p-2 text-right">${fmt(res.msBetween, 3)}</td>
            <td class="p-2 text-right font-bold">${fmt(res.fVal, 3)}</td>
            <td class="p-2 text-right font-bold">${fmtPRaw(res.pVal)}</td>
          </tr>
          <tr class="border-b border-slate-200">
            <td class="p-2 font-medium">Within Groups</td>
            <td class="p-2 text-right">${fmt(res.ssWithin, 3)}</td>
            <td class="p-2 text-right">${res.dfWithin}</td>
            <td class="p-2 text-right">${fmt(res.msWithin, 3)}</td>
            <td class="p-2 text-right">-</td>
            <td class="p-2 text-right">-</td>
          </tr>
          <tr class="bg-slate-50 font-bold">
            <td class="p-2">Total</td>
            <td class="p-2 text-right">${fmt(res.ssTotal, 3)}</td>
            <td class="p-2 text-right">${res.dfTotal}</td>
            <td class="p-2 text-right">-</td><td class="p-2 text-right">-</td><td class="p-2 text-right">-</td>
          </tr>
        </tbody>
      </table>
    `;

    appState.currentQuiz = [
      {
        q: "1. จากตาราง ANOVA องศาอิสระระหว่างกลุ่ม (df Between) และภายในกลุ่ม (df Within) มีค่าเท่าใด?",
        options: [
          `df Between = ${res.dfBetween}, df Within = ${res.dfWithin}`,
          `df Between = ${res.dfWithin}, df Within = ${res.dfTotal}`,
          `df Between = ${res.k}, df Within = ${res.N}`
        ],
        correct: 0,
        explanation: `df ระหว่างกลุ่ม = k - 1 = ${res.dfBetween} และ df ภายในกลุ่ม = N - k = ${res.dfWithin}`
      },
      {
        q: "2. ค่าสถิติ F ในตารางคำนวณมาจากข้อใด?",
        options: [
          `นำ MS Between (${fmt(res.msBetween, 3)}) หารด้วย MS Within (${fmt(res.msWithin, 3)}) ได้ F = ${fmt(res.fVal, 3)}`,
          `นำ SS Between (${fmt(res.ssBetween, 3)}) หารด้วย SS Total (${fmt(res.ssTotal, 3)})`
        ],
        correct: 0,
        explanation: `ค่า F คำนวณจาก MS ระหว่างกลุ่ม ÷ MS ภายในกลุ่ม = ${fmt(res.msBetween, 3)} ÷ ${fmt(res.msWithin, 3)} = ${fmt(res.fVal, 3)}`
      },
      {
        q: "3. จากค่า Sig. ในตาราง ข้อใดสรุปผลและระบุขั้นตอนต่อไปได้ถูกต้อง?",
        options: [
          "มีอย่างน้อย 1 คู่ที่แตกต่างกันอย่างมีนัยสำคัญ (.05) และต้องทดสอบรายคู่ (Post-hoc) ต่อ",
          "ทุกกลุ่มมีค่าเฉลี่ยไม่แตกต่างกัน จึงไม่ต้องทดสอบรายคู่ (Post-hoc)"
        ],
        correct: isSig ? 0 : 1,
        explanation: `ค่า Sig. = ${fmtPRaw(res.pVal)} ซึ่ง${isSig ? "น้อยกว่า 0.05 แสดงว่ามีอย่างน้อย 1 คู่ต่างกัน ต้องดูตาราง Post-hoc ต่อ" : "มากกว่าหรือเท่ากับ 0.05 แสดงว่าไม่แตกต่างกัน"}`
      }
    ];
  } else if (topic === "onesample") {
    const dv = appState.dvs[Math.floor(Math.random() * appState.dvs.length)];
    const mu0 = dv.defaultTestVal || 3.5;
    const vals = appState.data.map((r) => clampLikert(r[dv.id]));
    const res = StatMath.oneSampleTTest(vals, mu0, 0.05);
    const isSig = res.p < 0.05;
    const meaning = interpretMean5Level(res.mean);

    document.getElementById("quizScenarioTag").textContent = "โจทย์ฝึกอ่านตาราง: One-Sample t-test";
    document.getElementById("quizScenarioTitle").textContent = `ทดสอบค่าเฉลี่ย "${dv.label}" กับเกณฑ์มาตรฐาน Test Value = ${fmt(mu0, 2)} คะแนน`;

    tableContainer.innerHTML = `
      <table class="w-full text-left border-collapse text-xs bg-white border border-slate-300">
        <thead>
          <tr class="bg-slate-100 border-b border-slate-400">
            <th class="p-2">ตัวแปรตาม (5 ระดับ)</th><th class="p-2 text-right">n</th><th class="p-2 text-right"><span class="xbar">x</span> (Mean)</th><th class="p-2 text-right">S.D.</th><th class="p-2 text-right">t</th><th class="p-2 text-right">df</th><th class="p-2 text-right">Sig. (2-tailed)</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td class="p-2 font-medium">${dv.label}</td>
            <td class="p-2 text-right">${res.n}</td>
            <td class="p-2 text-right">${fmt(res.mean, 2)}</td>
            <td class="p-2 text-right">${fmt(res.sd, 2)}</td>
            <td class="p-2 text-right font-bold">${fmt(res.t, 3)}</td>
            <td class="p-2 text-right">${res.df}</td>
            <td class="p-2 text-right font-bold">${fmtPRaw(res.p)}</td>
          </tr>
        </tbody>
      </table>
    `;

    appState.currentQuiz = [
      {
        q: `1. ค่าเฉลี่ย (${fmt(res.mean, 2)}) แปลผลตามเกณฑ์ 5 ระดับได้ในระดับใด และต่างจากเกณฑ์ (${fmt(mu0, 2)}) อยู่เท่าใด?`,
        options: [
          `อยู่ในระดับ "${meaning.text}" และมีผลต่างค่าเฉลี่ย (Mean Difference) = ${fmt(res.diff, 2)} คะแนน`,
          `อยู่ในระดับ "น้อยที่สุด" และมีผลต่างค่าเฉลี่ย = ${fmt(res.sd, 2)} คะแนน`
        ],
        correct: 0,
        explanation: `ค่าเฉลี่ย ${fmt(res.mean, 2)} อยู่ในระดับ${meaning.text} และ Mean Difference = ${fmt(res.mean, 2)} − ${fmt(mu0, 2)} = ${fmt(res.diff, 2)}`
      },
      {
        q: "2. จากค่า Sig. (2-tailed) ข้อใดแปลผลที่ระดับนัยสำคัญ .05 ได้ถูกต้อง?",
        options: [
          `ค่าเฉลี่ยแตกต่างจากเกณฑ์ ${fmt(mu0, 2)} คะแนน อย่างมีนัยสำคัญทางสถิติที่ระดับ .05`,
          `ค่าเฉลี่ยไม่แตกต่างจากเกณฑ์ ${fmt(mu0, 2)} คะแนน ทางสถิติ`
        ],
        correct: isSig ? 0 : 1,
        explanation: `ค่า Sig. = ${fmtPRaw(res.p)} ซึ่ง${isSig ? "น้อยกว่า 0.05 จึงปฏิเสธ H₀" : "มากกว่าหรือเท่ากับ 0.05 จึงยอมรับ H₀"}`
      }
    ];
  } else {
    const iv = appState.ivs[0];
    const N = appState.data.length;
    const counts = iv.options.map((opt) => ({
      opt,
      c: appState.data.filter((r) => r[iv.id] === opt).length
    }));
    const sorted = [...counts].sort((a, b) => b.c - a.c);
    const top = sorted[0];
    const topPct = N > 0 ? (top.c / N) * 100 : 0;

    document.getElementById("quizScenarioTag").textContent = "โจทย์ฝึกอ่านตาราง: ความถี่และร้อยละ";
    document.getElementById("quizScenarioTitle").textContent = `ตารางแจกแจงความถี่และร้อยละ จำแนกตาม "${iv.label}"`;

    tableContainer.innerHTML = `
      <table class="w-full text-left border-collapse text-xs bg-white border border-slate-300">
        <thead>
          <tr class="bg-slate-100 border-b border-slate-400">
            <th class="p-2">${iv.label}</th><th class="p-2 text-right">ความถี่ (f)</th><th class="p-2 text-right">ร้อยละ (%)</th>
          </tr>
        </thead>
        <tbody>
          ${counts
            .map(
              (r) =>
                `<tr class="border-b border-slate-200"><td class="p-2">${r.opt}</td><td class="p-2 text-right">${r.c}</td><td class="p-2 text-right">${fmt((r.c / N) * 100, 2)}</td></tr>`
            )
            .join("")}
          <tr class="bg-slate-50 font-bold"><td class="p-2">รวม</td><td class="p-2 text-right">${N}</td><td class="p-2 text-right">100.00</td></tr>
        </tbody>
      </table>
    `;

    appState.currentQuiz = [
      {
        q: "1. กลุ่มที่มีจำนวนผู้ตอบมากที่สุดคือกลุ่มใด และคิดเป็นร้อยละเท่าใด?",
        options: [
          `กลุ่ม "${top.opt}" จำนวน ${top.c} คน (ร้อยละ ${fmt(topPct, 2)})`,
          `กลุ่ม "${sorted[sorted.length - 1].opt}" จำนวน ${N} คน (ร้อยละ 100.00)`
        ],
        correct: 0,
        explanation: `กลุ่ม "${top.opt}" มีความถี่สูงสุดคือ ${top.c} คน จากทั้งหมด ${N} คน คิดเป็นร้อยละ ${fmt(topPct, 2)}`
      }
    ];
  }

  qContainer.innerHTML = appState.currentQuiz
    .map(
      (item, qIdx) => `
      <div class="p-3.5 rounded-xl bg-slate-50 border border-slate-200 space-y-2">
        <p class="font-heading font-semibold text-xs sm:text-sm text-slate-800">${item.q}</p>
        <div class="space-y-1.5">
          ${item.options
            .map(
              (opt, oIdx) => `
              <label class="flex items-start gap-2 text-xs sm:text-sm text-slate-700 cursor-pointer p-1.5 rounded hover:bg-white transition">
                <input type="radio" name="quiz_q_${qIdx}" value="${oIdx}" class="mt-0.5 text-brand-600 focus:ring-brand-500" />
                <span>${opt}</span>
              </label>
            `
            )
            .join("")}
        </div>
      </div>
    `
    )
    .join("");
}

// ============================================================================
// 8. EVENT LISTENERS & INITIALIZATION
// ============================================================================

document.addEventListener("DOMContentLoaded", () => {
  const navTabs = document.querySelectorAll(".nav-tab");
  const tabPanels = document.querySelectorAll(".tab-panel");

  function switchMainTab(targetId) {
    navTabs.forEach((btn) => {
      if (btn.getAttribute("data-tab") === targetId) {
        btn.classList.add("active-tab");
      } else {
        btn.classList.remove("active-tab");
      }
    });
    tabPanels.forEach((panel) => {
      if (panel.id === targetId) {
        panel.classList.remove("hidden");
      } else {
        panel.classList.add("hidden");
      }
    });
  }

  navTabs.forEach((btn) => {
    btn.addEventListener("click", () => switchMainTab(btn.getAttribute("data-tab")));
  });

  document.getElementById("btnGoToCollect").addEventListener("click", () => switchMainTab("tab-collect"));
  document.getElementById("btnGoToAnalyze").addEventListener("click", () => switchMainTab("tab-analyze"));

  const statSubtabs = document.querySelectorAll(".stat-subtab");
  const statPanels = document.querySelectorAll(".stat-panel");

  statSubtabs.forEach((btn) => {
    btn.addEventListener("click", () => {
      const targetStat = btn.getAttribute("data-stat");
      statSubtabs.forEach((b) => b.classList.remove("active-stat-subtab"));
      btn.classList.add("active-stat-subtab");

      statPanels.forEach((p) => {
        if (p.id === `stat-${targetStat}`) {
          p.classList.remove("hidden");
        } else {
          p.classList.add("hidden");
        }
      });
    });
  });

  document.getElementById("toggleReadingGuide").addEventListener("change", (e) => {
    document.querySelectorAll(".reading-guide-box").forEach((box) => {
      if (e.target.checked) box.classList.remove("hidden");
      else box.classList.add("hidden");
    });
  });

  document.getElementById("btnLoadPreset").addEventListener("click", () => {
    const key = document.getElementById("presetSelect").value;
    if (confirm("ต้องการโหลดโจทย์ที่เลือกใหม่ (ข้อมูลที่แก้ไขค้างไว้จะถูกรีเซ็ตเป็นชุดนี้) ใช่หรือไม่?")) {
      loadPreset(key);
    }
  });
  document.getElementById("presetSelect").addEventListener("change", (e) => {
    loadPreset(e.target.value);
  });

  const btnBlank = document.getElementById("btnStartBlankProject");
  if (btnBlank) {
    btnBlank.addEventListener("click", () => {
      if (confirm("ต้องการเริ่มสร้างหัวข้อวิจัยใหม่ของตนเอง (ตารางข้อมูลเปล่า N=0) ใช่หรือไม่?")) {
        loadPreset("custom_blank");
      }
    });
  }

  // Add New Independent Variable (IV)
  document.getElementById("btnAddIv").addEventListener("click", () => {
    const labelInp = document.getElementById("newIvLabel");
    const optsInp = document.getElementById("newIvOptions");
    const label = labelInp.value.trim();
    const options = optsInp.value
      .split(",")
      .map((s) => s.trim())
      .filter(Boolean);

    if (!label || options.length < 2) {
      alert("กรุณาระบุชื่อตัวแปรต้น และตัวเลือกอย่างน้อย 2 กลุ่ม (คั่นด้วยเครื่องหมายจุลภาค ,)");
      return;
    }

    const id = "iv_" + Date.now();
    appState.ivs.push({ id, label, options });
    appState.data.forEach((r) => {
      r[id] = options[Math.floor(Math.random() * options.length)];
    });
    labelInp.value = "";
    optsInp.value = "";
    saveStateToLocal();
    renderAll();
  });

  // Add New Dependent Variable (5-Level DV)
  document.getElementById("btnAddDv").addEventListener("click", () => {
    const labelInp = document.getElementById("newDvLabel");
    const label = labelInp.value.trim();

    if (!label) {
      alert("กรุณาระบุชื่อประเด็นหรือข้อคำถามของตัวแปรตาม (5 ระดับ)");
      return;
    }

    const id = "dv_" + Date.now();
    appState.dvs.push({ id, label, unit: "ระดับ (1-5)", defaultTestVal: 3.5 });
    appState.data.forEach((r) => {
      r[id] = clampLikert(3 + Math.round((Math.random() - 0.3) * 2));
    });
    labelInp.value = "";
    saveStateToLocal();
    renderAll();
  });

  // Survey Form Submit (5-Level Radio Inputs)
  document.getElementById("surveyForm").addEventListener("submit", (e) => {
    e.preventDefault();
    const formData = new FormData(e.target);
    const newRow = {};
    appState.ivs.forEach((iv) => {
      newRow[iv.id] = formData.get(`iv_${iv.id}`);
    });
    appState.dvs.forEach((dv) => {
      newRow[dv.id] = clampLikert(formData.get(`dv_${dv.id}`));
    });
    appState.data.push(newRow);
    saveStateToLocal();
    renderProjectHeader();
    renderRawDataTable();
    runAllAnalyses();
  });

  document.getElementById("btnSimOne").addEventListener("click", () => {
    simulateRespondent();
    saveStateToLocal();
    renderProjectHeader();
    renderRawDataTable();
    runAllAnalyses();
  });

  document.getElementById("btnSimTen").addEventListener("click", () => {
    for (let i = 0; i < 10; i++) simulateRespondent();
    saveStateToLocal();
    renderProjectHeader();
    renderRawDataTable();
    runAllAnalyses();
  });

  document.getElementById("btnClearData").addEventListener("click", () => {
    if (confirm("ต้องการล้างข้อมูลทั้งหมดในตารางเพื่อเริ่มกรอกใหม่ใช่หรือไม่?")) {
      appState.data = [];
      saveStateToLocal();
      renderProjectHeader();
      renderRawDataTable();
      runAllAnalyses();
    }
  });

  // Import CSV from Google Forms / Excel
  const importInp = document.getElementById("inputImportCsv");
  if (importInp) {
    importInp.addEventListener("change", (e) => {
      const file = e.target.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = (evt) => {
        const text = evt.target.result;
        const lines = text
          .split(/\r?\n/)
          .map((l) => l.trim())
          .filter(Boolean);
        if (lines.length < 2) {
          alert("ไฟล์ CSV ต้องมีอย่างน้อย 1 แถวหัวตาราง และ 1 แถวข้อมูล");
          return;
        }
        const parseCsvLine = (line) =>
          line.split(",").map((cell) => cell.trim().replace(/^"|"$/g, ""));

        const headers = parseCsvLine(lines[0]);
        // Skip first column if it is ID/Timestamp
        const startCol =
          headers[0].toLowerCase() === "id" ||
          headers[0].includes("ประทับเวลา") ||
          headers[0].toLowerCase().includes("timestamp")
            ? 1
            : 0;

        const importedRows = [];
        for (let r = 1; r < lines.length; r++) {
          const cols = parseCsvLine(lines[r]);
          if (cols.length <= startCol) continue;
          const rowObj = {};
          appState.ivs.forEach((iv, idx) => {
            const rawVal = cols[startCol + idx] || iv.options[0];
            if (!iv.options.includes(rawVal) && rawVal) {
              iv.options.push(rawVal);
            }
            rowObj[iv.id] = rawVal || iv.options[0];
          });
          appState.dvs.forEach((dv, idx) => {
            const rawNum = cols[startCol + appState.ivs.length + idx];
            rowObj[dv.id] = clampLikert(rawNum);
          });
          importedRows.push(rowObj);
        }

        if (importedRows.length > 0) {
          appState.data = importedRows;
          saveStateToLocal();
          renderAll();
          alert(`นำเข้าข้อมูลสำเร็จ จำนวน ${importedRows.length} คน!`);
        }
        importInp.value = "";
      };
      reader.readAsText(file, "UTF-8");
    });
  }

  document.getElementById("btnExportCsv").addEventListener("click", () => {
    const headers = [
      "ID",
      ...appState.ivs.map((v) => v.label),
      ...appState.dvs.map((v) => `${v.label} (1-5)`)
    ];
    const rows = appState.data.map((r, idx) => [
      idx + 1,
      ...appState.ivs.map((v) => `"${r[v.id] ?? ""}"`),
      ...appState.dvs.map((v) => r[v.id] ?? 3)
    ]);
    const csvContent =
      "\uFEFF" + [headers.join(","), ...rows.map((e) => e.join(","))].join("\n");
    const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "statlab_5level_dataset.csv";
    a.click();
    URL.revokeObjectURL(url);
  });

  document.getElementById("freqVarSelect").addEventListener("change", renderFrequencyAnalysis);

  document.getElementById("oneSampleDvSelect").addEventListener("change", (e) => {
    const dv = appState.dvs.find((v) => v.id === e.target.value);
    if (dv && dv.defaultTestVal !== undefined) {
      document.getElementById("oneSampleTestVal").value = dv.defaultTestVal;
    }
    renderOneSampleAnalysis();
  });
  document.getElementById("oneSampleTestVal").addEventListener("input", renderOneSampleAnalysis);
  document.getElementById("oneSampleAlpha").addEventListener("change", renderOneSampleAnalysis);

  document.getElementById("indepIvSelect").addEventListener("change", () => {
    updateIndepGroupSelectors();
    renderIndepTTestAnalysis();
  });
  document.getElementById("indepGroup1Select").addEventListener("change", renderIndepTTestAnalysis);
  document.getElementById("indepGroup2Select").addEventListener("change", renderIndepTTestAnalysis);
  document.getElementById("indepDvSelect").addEventListener("change", renderIndepTTestAnalysis);

  document.getElementById("anovaIvSelect").addEventListener("change", renderAnovaAnalysis);
  document.getElementById("anovaDvSelect").addEventListener("change", renderAnovaAnalysis);

  document.getElementById("quizTopicSelect").addEventListener("change", generateQuiz);
  document.getElementById("btnNewQuiz").addEventListener("click", generateQuiz);

  document.getElementById("btnCheckQuiz").addEventListener("click", () => {
    if (!appState.currentQuiz) return;
    let score = 0;
    const details = [];

    appState.currentQuiz.forEach((item, idx) => {
      const selected = document.querySelector(`input[name="quiz_q_${idx}"]:checked`);
      const userAns = selected ? Number(selected.value) : -1;
      const isCorrect = userAns === item.correct;
      if (isCorrect) score++;
      details.push(`
        <div class="p-2.5 rounded-lg border ${
          isCorrect ? "bg-emerald-50 border-emerald-200 text-emerald-950" : "bg-rose-50 border-rose-200 text-rose-950"
        }">
          <div class="font-bold">${isCorrect ? "✅ ข้อที่ " + (idx + 1) + " ถูกต้อง!" : "❌ ข้อที่ " + (idx + 1) + " ยังไม่ถูกต้อง"}</div>
          <div class="mt-0.5 text-slate-700"><strong>คำอธิบาย:</strong> ${item.explanation}</div>
        </div>
      `);
    });

    const fb = document.getElementById("quizFeedbackBox");
    fb.className = "rounded-xl p-4 text-xs sm:text-sm space-y-2 bg-slate-100 border border-slate-300 block";
    fb.innerHTML = `
      <div class="font-heading font-bold text-base text-slate-800">
        คะแนนของคุณ: ${score} / ${appState.currentQuiz.length} ข้อ
      </div>
      <div class="space-y-2">${details.join("")}</div>
    `;
  });

  const modal = document.getElementById("titleModal");
  document.getElementById("btnEditProjectTitle").addEventListener("click", () => {
    document.getElementById("inputProjectTitle").value = appState.title;
    document.getElementById("inputProjectDesc").value = appState.description;
    modal.classList.remove("hidden");
    modal.classList.add("flex");
  });
  document.getElementById("btnCancelTitleModal").addEventListener("click", () => {
    modal.classList.add("hidden");
    modal.classList.remove("flex");
  });
  document.getElementById("btnSaveTitleModal").addEventListener("click", () => {
    appState.title = document.getElementById("inputProjectTitle").value.trim() || appState.title;
    appState.description = document.getElementById("inputProjectDesc").value.trim() || appState.description;
    saveStateToLocal();
    renderProjectHeader();
    modal.classList.add("hidden");
    modal.classList.remove("flex");
  });

  // Load from localStorage first; if not found, load default preset
  if (loadStateFromLocal()) {
    renderAll();
  } else {
    loadPreset("student_cafe");
  }
});
