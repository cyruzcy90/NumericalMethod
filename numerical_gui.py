import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def preprocess_equation(expr):
    import re
    
    expr = expr.replace(" ", "")
    
    functions = ["sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", 
                 "exp", "log", "log10", "sqrt", "abs", "floor", "ceil", "fabs", "arcsin", 
                 "arccos", "arctan", "ln", "log2", "deg", "rad", "trunc"]

    # Handle implicit multiplication:
    # 1. Number followed by 'x' or function
    expr = re.sub(r'(\d)([x(])', r'\1*\2', expr)
    # 2. 'x' or ')' followed by number or 'x' or function
    expr = re.sub(r'([x\)])(\d|[x(])', r'\1*\2', expr)
    # 3. Number or 'x' followed by function name
    for func in functions:
        expr = re.sub(rf'(\d|x)({func}\()', r'\1*\2', expr)
    # 4. 'x' followed by '(', e.g., x(x+1)
    expr = re.sub(r'(x)(\()', r'\1*\2', expr)
    # 5. ')' followed by '(', e.g., (x+1)(x-1)
    expr = re.sub(r'(\))(\()', r'\1*\2', expr)
    # 6. 'x' followed by 'x'
    expr = re.sub(r'(x)(x)', r'\1*\2', expr)

    return expr

def safe_eval(expr, x):
 
    try:
        expr = preprocess_equation(expr)

        namespace = {
            "x": x,
            # Trigonometric functions
            "sin": math.sin, "cos": math.cos, "tan": math.tan,
            "asin": math.asin, "acos": math.acos, "atan": math.atan,
            "arcsin": math.asin, "arccos": math.acos, "arctan": math.atan,
            # Hyperbolic functions
            "sinh": math.sinh, "cosh": math.cosh, "tanh": math.tanh,
            # Exponential and logarithmic
            "exp": math.exp, "log": math.log, "ln": math.log, 
            "log10": math.log10, "log2": math.log2,
            # Root and power
            "sqrt": math.sqrt,
            # Absolute and rounding
            "abs": abs, "fabs": math.fabs,
            "floor": math.floor, "ceil": math.ceil,
            "trunc": math.trunc,
            # Constants
            "pi": math.pi, "e": math.e,
            # Additional utilities
            "deg": math.degrees, "rad": math.radians,
        }

     
        result = eval(expr, {"__builtins__": {}}, namespace)
        return result
    except:
        return None

# ─────────────────────────────────────────────
#  NUMERICAL METHODS
# ─────────────────────────────────────────────
def incremental_method(f, x0, dx=0.1, steps=100):
    rows, root = [], None
    xa = x0
    for _ in range(steps):
        xb = xa + dx
        fa, fb = safe_eval(f, xa), safe_eval(f, xb)
        if fa is None or fb is None:
            break
        rows.append((round(xa,6), round(xb,6), round(fa,6), round(fb,6), round(fa*fb,6)))
        if fa * fb < 0:
            root = (xa + xb) / 2
            break
        xa = xb
    return rows, root

def bisection_method(f, a, b, tol=1e-6, maxiter=100):
    rows, root = [], None
    fa, fb = safe_eval(f, a), safe_eval(f, b)
    if fa * fb > 0:
        return rows, None
    for i in range(1, maxiter+1):
        c = (a + b) / 2 
        fc = safe_eval(f, c)
        rows.append((i, round(a,6), round(b,6), round(c,6), round(fc,6)))
        if abs(fc) < tol or (b - a) / 2 < tol:
            root = c; break
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return rows, root

def regula_falsi(f, a, b, tol=1e-6, maxiter=100):
    rows, root = [], None
    fa, fb = safe_eval(f, a), safe_eval(f, b)
    if fa * fb > 0:
        return rows, None
    for i in range(1, maxiter+1):
        if abs(fb - fa) < 1e-15:
            break
        c = b - fb * (b - a) / (fb - fa)
        fc = safe_eval(f, c)
        rows.append((i, round(a,6), round(b,6), round(c,6), round(fc,6)))
        if abs(fc) < tol:
            root = c; break
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return rows, root

def newton_raphson(f, x0, tol=1e-6, maxiter=100):
    rows, root, h = [], None, 1e-7
    x = x0
    for i in range(1, maxiter+1):
        fx = safe_eval(f, x)
        dfx = (safe_eval(f, x+h) - safe_eval(f, x-h)) / (2*h)
        if dfx is None or abs(dfx) < 1e-15:
            break
        x1 = x - fx / dfx
        rows.append((i, round(x,6), round(fx,6), round(dfx,6), round(x1,6)))
        if abs(x1 - x) < tol:
            root = x1; break
        x = x1
    return rows, root

def secant_method(f, x0, x1, tol=1e-6, maxiter=100):
    rows, root = [], None
    for i in range(1, maxiter+1):
        f0, f1 = safe_eval(f, x0), safe_eval(f, x1)
        if abs(f1 - f0) < 1e-15:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        rows.append((i, round(x0,6), round(x1,6), round(x2,6), round(safe_eval(f,x2),6)))
        if abs(x2 - x1) < tol:
            root = x2; break
        x0, x1 = x1, x2
    return rows, root

# ─────────────────────────────────────────────
#  MATRIX OPERATIONS 
# ─────────────────────────────────────────────
def mat_add(A, B):
    n, m = len(A), len(A[0])
    return [[A[i][j]+B[i][j] for j in range(m)] for i in range(n)]

def mat_mul(A, B):
    n, m, p = len(A), len(A[0]), len(B[0])
    return [[sum(A[i][k]*B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]

def mat_transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]

def mat_det(A):
    n = len(A)
    if n == 1: return A[0][0]
    if n == 2: return A[0][0]*A[1][1] - A[0][1]*A[1][0]
    d = 0
    for c in range(n):
        minor = [[A[r][cc] for cc in range(n) if cc!=c] for r in range(1,n)]
        d += ((-1)**c) * A[0][c] * mat_det(minor)
    return d

def mat_cofactor(A):
    n = len(A)
    C = []
    for i in range(n):
        row = []
        for j in range(n):
            minor = [[A[r][c] for c in range(n) if c!=j] for r in range(n) if r!=i]
            row.append(((-1)**(i+j)) * mat_det(minor))
        C.append(row)
    return C

def mat_adjoint(A):
    return mat_transpose(mat_cofactor(A))

def mat_inverse(A):
    d = mat_det(A)
    if abs(d) < 1e-12:
        return None
    adj = mat_adjoint(A)
    n = len(A)
    return [[adj[i][j]/d for j in range(n)] for i in range(n)]

def mat_power(A, p):
    n = len(A)
    R = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    for _ in range(p):
        R = mat_mul(R, A)
    return R

def solve_linear(A, b):
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(M[r][col]))
        M[col], M[pivot] = M[pivot], M[col]
        if abs(M[col][col]) < 1e-12:
            return None
        for row in range(n):
            if row != col:
                f = M[row][col] / M[col][col]
                M[row] = [M[row][k] - f*M[col][k] for k in range(n+1)]
    return [M[i][n] / M[i][i] for i in range(n)]

def parse_matrix(text):
    try:
        rows = [r.strip() for r in text.strip().split('\n') if r.strip()]
        return [[float(v) for v in r.replace(',', ' ').split()] for r in rows]
    except:
        return None

def fmt_matrix(M, decimals=4):
    if M is None: return "No result (singular or invalid)"
    lines = []
    for row in M:
        lines.append("  " + "  ".join(f"{v:>{10}.{decimals}f}" for v in row))
    return "\n".join(lines)

# ─────────────────────────────────────────────
#  MAIN APPLICATION
# ─────────────────────────────────────────────
class NumericalApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Numerical Methods Solver")
        self.geometry("1100x720")
        self.configure(bg="#0f1117")
        self._setup_styles()
        self._build_ui()

    def _setup_styles(self):
        self.style = ttk.Style(self)
        self.style.theme_use("clam")
        BG, FG, ACC = "#0f1117", "#e2e8f0", "#38bdf8"
        CARD = "#1e2330"
        self.style.configure("TNotebook", background=BG, borderwidth=0, tabmargins=[0,0,0,0])
        self.style.configure("TNotebook.Tab", background="#1a1f2e", foreground="#94a3b8",
                             padding=[18,8], font=("Courier New",10,"bold"), borderwidth=0)
        self.style.map("TNotebook.Tab",
                       background=[("selected", ACC)],
                       foreground=[("selected", "#0f1117")])
        self.style.configure("TFrame", background=BG)
        self.style.configure("Card.TFrame", background=CARD)
        self.style.configure("TLabel", background=BG, foreground=FG, font=("Courier New",10))
        self.style.configure("Title.TLabel", background=BG, foreground=ACC,
                             font=("Courier New",13,"bold"))
        self.style.configure("TButton", background=ACC, foreground="#0f1117",
                             font=("Courier New",10,"bold"), borderwidth=0, padding=[10,5])
        self.style.map("TButton", background=[("active","#0ea5e9")])
        self.style.configure("TCombobox", fieldbackground=CARD, background=CARD,
                             foreground=FG, font=("Courier New",10))

    # ── TOP HEADER ──────────────────────────────
    def _build_ui(self):
        hdr = tk.Frame(self, bg="#0f1117", pady=10)
        hdr.pack(fill="x", padx=20)
        tk.Label(hdr, text="NUMERICAL METHODS SOLVER", bg="#0f1117",
                 fg="#38bdf8", font=("Courier New",16,"bold")).pack(side="left")
        tk.Label(hdr, text="v1.0  |  Root Finding + Matrix Operations",
                 bg="#0f1117", fg="#475569", font=("Courier New",9)).pack(side="left", padx=12)

        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=10, pady=(0,10))

        self._build_root_tab(nb)
        self._build_matrix_tab(nb)

    # ── ROOT FINDING TAB ────────────────────────
    def _build_root_tab(self, nb):
        tab = ttk.Frame(nb)
        nb.add(tab, text="  Root Finding  ")

        left = tk.Frame(tab, bg="#1e2330", padx=14, pady=14)
        left.pack(side="left", fill="y", padx=(8,4), pady=8)

        def lbl(parent, text, bold=False):
            f = ("Courier New",10,"bold") if bold else ("Courier New",10)
            tk.Label(parent, text=text, bg="#1e2330", fg="#94a3b8", font=f,
                     anchor="w").pack(fill="x", pady=(6,1))

        lbl(left, "f(x) equation:", bold=True)
        self.fx_entry = tk.Entry(left, bg="#0f1117", fg="#f8fafc", font=("Courier New",11),
                                  insertbackground="#38bdf8", width=30, relief="flat", bd=4)
        self.fx_entry.insert(0, "x**4 - 5*x + 4*x**2")
        self.fx_entry.pack(fill="x")

        # Help button
        help_frame = tk.Frame(left, bg="#1e2330")
        help_frame.pack(fill="x", pady=(2,0))
        tk.Label(help_frame, text="Examples:", bg="#1e2330", fg="#64748b", 
                 font=("Courier New",8, "bold")).pack(side="left")
        tk.Label(help_frame, text="x**2-4, sin(x)-x, x**3-1, exp(x)-3, sqrt(x)-2",
                 bg="#1e2330", fg="#64748b", font=("Courier New",8), wraplength=250, justify="left").pack(side="left", padx=(2,0))

        # Supported functions info
        info_text = """Supported: +,-,*,/,**,x,()
sin,cos,tan,asin,acos,atan
sinh,cosh,tanh
exp,log,ln,log10,log2,sqrt
abs,floor,ceil,pi,e"""

        info_label = tk.Label(left, text=info_text, bg="#1e2330", fg="#475569", 
                              font=("Courier New",7), justify="left")
        info_label.pack(fill="x", pady=(6,0))

        lbl(left, "Method:")
        self.method_var = tk.StringVar(value="Auto Detect")
        methods = ["Auto Detect","Incremental","Bisection","Regula-Falsi","Newton-Raphson","Secant"]
        self.method_cb = ttk.Combobox(left, values=methods, textvariable=self.method_var,
                                       state="readonly", width=28)
        self.method_cb.pack(fill="x")

        lbl(left, "x0 / a:")
        self.x0e = tk.Entry(left, bg="#0f1117", fg="#f8fafc", font=("Courier New",11),
                             insertbackground="#38bdf8", relief="flat", bd=4, width=30)
        self.x0e.insert(0, "1"); self.x0e.pack(fill="x")

        lbl(left, "x1 / b (for bracketing):")
        self.x1e = tk.Entry(left, bg="#0f1117", fg="#f8fafc", font=("Courier New",11),
                             insertbackground="#38bdf8", relief="flat", bd=4, width=30)
        self.x1e.insert(0, "2"); self.x1e.pack(fill="x")

        lbl(left, "Tolerance:")
        self.tole = tk.Entry(left, bg="#0f1117", fg="#f8fafc", font=("Courier New",11),
                              insertbackground="#38bdf8", relief="flat", bd=4, width=30)
        self.tole.insert(0, "1e-6"); self.tole.pack(fill="x")

        lbl(left, "Step dx (Incremental):")
        self.dxe = tk.Entry(left, bg="#0f1117", fg="#f8fafc", font=("Courier New",11),
                             insertbackground="#38bdf8", relief="flat", bd=4, width=30)
        self.dxe.insert(0, "0.1"); self.dxe.pack(fill="x")

        ttk.Button(left, text="SOLVE", command=self._solve_root).pack(fill="x", pady=(14,4))
        ttk.Button(left, text="CLEAR", command=lambda: [
            self.root_table.delete(*self.root_table.get_children()),
            self.root_result.config(text="")
        ]).pack(fill="x")

        right = tk.Frame(tab, bg="#0f1117")
        right.pack(side="right", fill="both", expand=True, padx=(4,8), pady=8)

        self.root_result = tk.Label(right, text="", bg="#0f1117", fg="#4ade80",
                                     font=("Courier New",12,"bold"))
        self.root_result.pack(anchor="w", padx=4, pady=(4,2))

        # table frame
        tf = tk.Frame(right, bg="#0f1117")
        tf.pack(fill="x", padx=4)

        self.root_tree_frame = tf
        self.root_table = None

        # graph
        self.fig, self.ax = plt.subplots(figsize=(5,2.8), facecolor="#0f1117")
        self.ax.set_facecolor("#1e2330")
        self.canvas_plot = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas_plot.get_tk_widget().pack(fill="both", expand=True, padx=4, pady=4)

    def _make_table(self, parent, cols):
        for w in parent.winfo_children():
            w.destroy()
        s = ttk.Style()
        s.configure("Root.Treeview", background="#1e2330", foreground="#e2e8f0",
                    fieldbackground="#1e2330", font=("Courier New",9), rowheight=22)
        s.configure("Root.Treeview.Heading", background="#0f1117", foreground="#38bdf8",
                    font=("Courier New",9,"bold"))
        s.map("Root.Treeview", background=[("selected","#38bdf8")], foreground=[("selected","#0f1117")])
        tree = ttk.Treeview(parent, columns=cols, show="headings", height=8, style="Root.Treeview")
        for c in cols:
            tree.heading(c, text=c)
            tree.column(c, width=max(80, int(560/len(cols))), anchor="center")
        sb = ttk.Scrollbar(parent, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=sb.set)
        tree.pack(side="left", fill="x", expand=True)
        sb.pack(side="right", fill="y")
        return tree

    def _solve_root(self):
        f  = self.fx_entry.get().strip()
        m  = self.method_var.get()
        try:
            x0 = float(self.x0e.get())
            x1 = float(self.x1e.get())
            tol = float(self.tole.get())
            dx  = float(self.dxe.get())
        except:
            messagebox.showerror("Input Error", "Check numeric inputs."); return

        if not f:
            messagebox.showerror("Error", "Please enter an equation."); return

        # Test equation with x0 to catch syntax errors early
        fa = safe_eval(f, x0)
        if fa is None:
            messagebox.showerror("Error", 
                f"Invalid equation. Check syntax.\n"
                f"Supported: +,-,*,/,**, x, ()\nsin, cos, tan, exp, log, sqrt, abs, etc.")
            return

        fb = safe_eval(f, x1)
        if fb is None:
            messagebox.showerror("Error", "Cannot evaluate equation at x1. Check equation syntax."); return

        bracket = fa is not None and fb is not None and fa*fb < 0

        if m == "Auto Detect":
            if bracket:  m = "Bisection"
            else:        m = "Newton-Raphson"

        rows, root, cols = [], None, []

        if m == "Incremental":
            cols = ["xa","xb","f(xa)","f(xb)","f(xa)*f(xb)"]
            rows, root = incremental_method(f, x0, dx)
        elif m == "Bisection":
            cols = ["Iter","a","b","c","f(c)"]
            rows, root = bisection_method(f, x0, x1, tol)
        elif m == "Regula-Falsi":
            cols = ["Iter","a","b","c","f(c)"]
            rows, root = regula_falsi(f, x0, x1, tol)
        elif m == "Newton-Raphson":
            cols = ["Iter","x","f(x)","f'(x)","x_new"]
            rows, root = newton_raphson(f, x0, tol)
        elif m == "Secant":
            cols = ["Iter","x0","x1","x2","f(x2)"]
            rows, root = secant_method(f, x0, x1, tol)

        self.root_table = self._make_table(self.root_tree_frame, cols)
        for r in rows:
            self.root_table.insert("", "end", values=r)

        if root is not None:
            self.root_result.config(text=f"  Root found [{m}]:  x = {root:.8f}   f(x) = {safe_eval(f,root):.2e}", fg="#4ade80")
        else:
            self.root_result.config(text=f"  [{m}] No root found in given interval.", fg="#f87171")

        self._plot_function(f, x0, x1, root, m)

    def _plot_function(self, f, a, b, root, method):
        self.ax.clear()
        self.ax.set_facecolor("#1e2330")
        span = max(abs(b-a)*1.5, 2)
        mid  = (a+b)/2
        xs   = [mid - span + span*2*i/400 for i in range(401)]
        ys   = [safe_eval(f,x) for x in xs]
        valid = [(x,y) for x,y in zip(xs,ys) if y is not None and abs(y)<1e6]
        if valid:
            vx, vy = zip(*valid)
            self.ax.plot(vx, vy, color="#38bdf8", linewidth=2, label=f"f(x)")
        self.ax.axhline(0, color="#475569", linewidth=1)
        self.ax.axvline(0, color="#475569", linewidth=1)
        if root is not None:
            self.ax.plot(root, safe_eval(f,root), "o", color="#4ade80", markersize=9,
                         zorder=5, label=f"Root x={root:.5f}")
        self.ax.legend(facecolor="#1e2330", edgecolor="#38bdf8", labelcolor="#e2e8f0",
                       fontsize=8)
        self.ax.set_title(f"{method}  |  f(x) = {f}", color="#94a3b8", fontsize=9)
        self.ax.tick_params(colors="#475569")
        for sp in self.ax.spines.values(): sp.set_edgecolor("#2d3748")
        self.fig.tight_layout()
        self.canvas_plot.draw()

    # ── MATRIX TAB ──────────────────────────────
    def _build_matrix_tab(self, nb):
        tab = ttk.Frame(nb)
        nb.add(tab, text="  Matrix Operations  ")

        left = tk.Frame(tab, bg="#1e2330", padx=14, pady=14, width=320)
        left.pack(side="left", fill="y", padx=(8,4), pady=8)
        left.pack_propagate(False)

        def lbl(txt, bold=False):
            f = ("Courier New",10,"bold") if bold else ("Courier New",10)
            tk.Label(left, text=txt, bg="#1e2330", fg="#94a3b8", font=f, anchor="w").pack(fill="x", pady=(6,1))

        lbl("Operation:", bold=True)
        ops = ["Add","Multiply","Transpose","Determinant","Adjoint","Inverse","Power","Solve Ax=b"]
        self.mat_op = tk.StringVar(value="Add")
        ttk.Combobox(left, values=ops, textvariable=self.mat_op, state="readonly", width=28).pack(fill="x")

        lbl("Matrix A  (rows, space-separated):", bold=True)
        self.mat_a = tk.Text(left, height=5, bg="#0f1117", fg="#f8fafc",
                              font=("Courier New",10), insertbackground="#38bdf8",
                              relief="flat", bd=4, width=30)
        self.mat_a.insert("1.0","1 2 3\n4 5 6\n7 8 9")
        self.mat_a.pack(fill="x")

        lbl("Matrix B  (for Add / Multiply / Solve b-vector):", bold=True)
        self.mat_b = tk.Text(left, height=5, bg="#0f1117", fg="#f8fafc",
                              font=("Courier New",10), insertbackground="#38bdf8",
                              relief="flat", bd=4, width=30)
        self.mat_b.insert("1.0","9 8 7\n6 5 4\n3 2 1")
        self.mat_b.pack(fill="x")

        lbl("Power (n)  for Matrix^n:")
        self.pow_e = tk.Entry(left, bg="#0f1117", fg="#f8fafc", font=("Courier New",11),
                               insertbackground="#38bdf8", relief="flat", bd=4, width=30)
        self.pow_e.insert(0, "2"); self.pow_e.pack(fill="x")

        ttk.Button(left, text="COMPUTE", command=self._compute_matrix).pack(fill="x", pady=(14,4))
        ttk.Button(left, text="CLEAR", command=lambda: self.mat_out.delete("1.0","end")).pack(fill="x")

        right = tk.Frame(tab, bg="#0f1117")
        right.pack(side="right", fill="both", expand=True, padx=(4,8), pady=8)

        tk.Label(right, text="Result", bg="#0f1117", fg="#38bdf8",
                 font=("Courier New",13,"bold")).pack(anchor="w", padx=6, pady=(4,2))

        self.mat_out = scrolledtext.ScrolledText(
            right, bg="#1e2330", fg="#f8fafc", font=("Courier New",11),
            relief="flat", bd=6, wrap="none", insertbackground="#38bdf8"
        )
        self.mat_out.pack(fill="both", expand=True, padx=6, pady=4)

    def _compute_matrix(self):
        op = self.mat_op.get()
        A  = parse_matrix(self.mat_a.get("1.0","end"))
        B  = parse_matrix(self.mat_b.get("1.0","end"))
        out = ""

        try:
            if op == "Add":
                if A is None or B is None: raise ValueError("Need A and B")
                R = mat_add(A, B)
                out = f"A + B =\n{fmt_matrix(R)}"

            elif op == "Multiply":
                if A is None or B is None: raise ValueError("Need A and B")
                R = mat_mul(A, B)
                out = f"A x B =\n{fmt_matrix(R)}"

            elif op == "Transpose":
                if A is None: raise ValueError("Need A")
                R = mat_transpose(A)
                out = f"Transpose(A) =\n{fmt_matrix(R)}"

            elif op == "Determinant":
                if A is None: raise ValueError("Need A")
                d = mat_det(A)
                out = f"det(A) = {d:.6f}"

            elif op == "Adjoint":
                if A is None: raise ValueError("Need A")
                R = mat_adjoint(A)
                out = f"Adjoint(A) =\n{fmt_matrix(R)}"

            elif op == "Inverse":
                if A is None: raise ValueError("Need A")
                R = mat_inverse(A)
                if R is None:
                    out = "Matrix is singular — no inverse."
                else:
                    out = f"Inverse(A) =\n{fmt_matrix(R)}"

            elif op == "Power":
                if A is None: raise ValueError("Need A")
                p = int(self.pow_e.get())
                R = mat_power(A, p)
                out = f"A^{p} =\n{fmt_matrix(R)}"

            elif op == "Solve Ax=b":
                if A is None or B is None: raise ValueError("Need A and b vector")
                b_vec = [row[0] for row in B]
                x = solve_linear(A, b_vec)
                if x is None:
                    out = "System has no unique solution (singular)."
                else:
                    out = "Solution x:\n" + "\n".join(
                        f"  x{i+1} = {v:.8f}" for i,v in enumerate(x))

        except Exception as e:
            out = f"Error: {e}"

        self.mat_out.delete("1.0","end")
        self.mat_out.insert("1.0", out)

# ─────────────────────────────────────────────
if __name__ == "__main__":
    app = NumericalApp()
