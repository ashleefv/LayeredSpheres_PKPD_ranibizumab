import os

import numpy as np
import pandas as pd

import itertools

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

import torch

import gpytorch
from gpytorch.kernels import MaternKernel, ScaleKernel, RBFKernel, LinearKernel, PeriodicKernel
from gpytorch.constraints.constraints import Interval
from gpytorch.likelihoods.gaussian_likelihood import GaussianLikelihood
from gpytorch.constraints import GreaterThan
from gpytorch.mlls import ExactMarginalLogLikelihood
from gpytorch.priors import GammaPrior

import botorch
from botorch.acquisition import ExpectedImprovement, UpperConfidenceBound, PosteriorMean
from botorch.acquisition.monte_carlo import SampleReducingMCAcquisitionFunction
from botorch.acquisition.monte_carlo import qExpectedImprovement
from botorch.acquisition.objective import ConstrainedMCObjective
from botorch.optim import optimize_acqf, optimize_acqf_mixed
import botorch.sampling.get_sampler
from botorch.sampling import SobolQMCNormalSampler
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP, MixedSingleTaskGP
from botorch.models.transforms.outcome import Standardize
from botorch.models.transforms.input import Normalize
from botorch.utils.sampling import draw_sobol_samples

import warnings
warnings.filterwarnings('ignore')

# plt.rcParams.update({'font.size':14})
torch.manual_seed(0)

def save_plot(fig, plot_name, df, condn_col, dpi = 600):
    """
    Saves a Matplotlib figure as a PNG under ../outputs/iteration_{n}/
    where n is the next iteration number (max existing + 1).
    The filename includes the plot name and the dose.
    """
    condn_val = float(df[condn_col].iloc[0])
    condn_str = f"{condn_val}mg"
    
    # Next iteration number
    next_iter = int(df["iteration"].max()) + 1
    
    # Path for outputs rel to helpers.py location
    script_dir = os.path.dirname(os.path.abspath(__file__))
#     output_dir = os.path.normpath(os.path.join(script_dir, "..", "outputs", f"iteration_{next_iter}"))
    output_dir = os.path.normpath(os.path.join(script_dir, "bo_outputs", f"iteration_{next_iter}"))
    os.makedirs(output_dir, exist_ok = True)

    filename = f"{plot_name.replace(' ', '_').lower()}_{condn_str}.png"
    save_path = os.path.join(output_dir, filename)

    fig.savefig(save_path, dpi = dpi, bbox_inches = "tight", pad_inches = 0.15)
    print(f"\nSaved: {save_path}")


def save_iteration_report(next_point, acq_value, pred_mean, std, target_col, df, condn_col, figs, ls_info, noise, std_ranges):
    
    next_iter = int(df["iteration"].max()) + 1
    script_dir = os.path.dirname(os.path.abspath(__file__))
#     output_dir = os.path.normpath(os.path.join(script_dir, "..", "outputs", f"iteration_{next_iter}"))
    output_dir = os.path.normpath(os.path.join(script_dir, "bo_outputs", f"iteration_{next_iter}"))

    os.makedirs(output_dir, exist_ok = True)
    
    pdf_path = os.path.join(output_dir, f"report_iter{next_iter}.pdf")
    
    with PdfPages(pdf_path) as pdf:
        
        # --- Text summary page ---
        fig_text, ax = plt.subplots(figsize = (8.5, 11))
        ax.axis('off')
        
        text = f"Iteration {next_iter}\n"
        text += "="*40 + "\n\n"

        text += "GP Diagnostics:\n"
        for k, v in ls_info.items():
            text += f"  {k} {v}\n"
        text += f"  Noise: {noise:.6g}\n\n"

        text += "Suggested next input:\n"
        for k, v in next_point.items():
            text += f"  {k} = {v:.6g}\n" if isinstance(v, float) else f"  {k} = {v}\n"
        text += f"\nEI at this point = {acq_value:.6g}\n\n"
        
        text += f"GP-predicted {target_col} = {pred_mean:.6g}\n"
        text += f"GP uncertainty (std) = {std:.6g}\n\n"
        
        # STD ranges per slice
        text += "STD ranges on slices:\n"
        for feature, std_min, std_max in std_ranges:
            text += f"  {feature}: [{std_min:.6g}  {std_max:.6g}]\n"
        text += "\n"        
        
        ax.text(0.1, 0.9, text, transform = ax.transAxes,
                fontsize = 12, verticalalignment = 'top', fontfamily = 'monospace')
        pdf.savefig(fig_text, bbox_inches = 'tight')
        plt.close(fig_text)
        
        for fig in figs:
            pdf.savefig(fig, bbox_inches = 'tight')
    
    print(f"\nReport saved: {pdf_path}")    
    
def get_highlight_points(df, target_col):
    """
    Returns:
        observed_color: str
        highlight_points: list of (df_high, color)
    """
    observed_color = "red"

    color_list = [
        "cyan", "magenta", "lime", "orange", "blue",
        "brown", "silver", "gold", "green", "chocolate", "orchid"
    ]

    highlight_points = []
    iteration_max = df["iteration"].max()

    if iteration_max == 0:
        idx = df[target_col].idxmax()
        df_high = df.loc[[idx]]
        highlight_points.append((df_high, "cyan"))
    else:
        for ds_id in range(0, iteration_max + 1):
            df_sub = df[df["iteration"] == ds_id]
            if len(df_sub) == 0:
                continue

            idx = df_sub[target_col].idxmax()
            df_high = df.loc[[idx]]
            color = color_list[ds_id % len(color_list)]
            highlight_points.append((df_high, color))

    return observed_color, highlight_points


def align_yzeros_pixel_exact(ax_left, ax_right, y_left=0.0, y_right=0.0,
                             ensure_zero_in_left=True, post_pad_right_max=None):
    """
    Align y=y_left on ax_left to y=y_right on ax_right at the exact same pixel height.
    Call this AFTER you finish plotting (including legends/tight_layout that may autoscale).

    Parameters
    ----------
    ax_left, ax_right : Matplotlib Axes
        ax_right should be created via ax_left.twinx().
    y_left, y_right : float
        Values to align (usually both 0.0).
    ensure_zero_in_left : bool
        Expand left y-limits a bit to include y_left if needed (recommended).
    post_pad_right_max : float or None
        If provided, ensures the right-axis top is >= this value to avoid clipping EI.
    """
    fig = ax_left.figure
    fig.canvas.draw()  # finalize transforms and autoscaled limits

    # Ensure y_left is inside left limits (alignment is ambiguous otherwise)
    y1_min, y1_max = ax_left.get_ylim()
    if ensure_zero_in_left and not (y1_min <= y_left <= y1_max):
        pad = 0.05 * (y1_max - y1_min + 1e-12)
        ax_left.set_ylim(min(y1_min, y_left) - pad, max(y1_max, y_left) + pad)
        fig.canvas.draw()
        y1_min, y1_max = ax_left.get_ylim()

    # Pixel row for y_left on the left axis
    ypix_target = ax_left.transData.transform((0.0, y_left))[1]

    # Left axis pixel span and fraction for y_left
    y1pix_min = ax_left.transData.transform((0.0, y1_min))[1]
    y1pix_max = ax_left.transData.transform((0.0, y1_max))[1]
    frac_left = (ypix_target - y1pix_min) / (y1pix_max - y1pix_min + 1e-12)

    # Keep the right data-span; shift its limits so y_right lands at the same fraction
    y2_min, y2_max = ax_right.get_ylim()
    span2 = (y2_max - y2_min) if (y2_max > y2_min) else 1.0
    new_y2_min = y_right - frac_left * span2
    new_y2_max = new_y2_min + span2
    ax_right.set_ylim(new_y2_min, new_y2_max)

    # Optional: ensure EI peak isn’t clipped
    if post_pad_right_max is not None:
        y2_min2, y2_max2 = ax_right.get_ylim()
        if y2_max2 < post_pad_right_max:
            ax_right.set_ylim(y2_min2, post_pad_right_max)

    fig.canvas.draw()
    
    

# Helper function that takes input data, kernel+likelihood specifications--> return trained model
def get_trained_GP(X, Y, kernel_type, cat_dims = None, noise_free = False, ls_floor = 0.20, ls_ceiling = 3.0, noise_prior_mean = 0.08):
    """
    This function is used to train a GP model based on the type of kernel that you select.
    This function will be used in modules 1,2,and 3. Save this function for later
    If you delete this cell accidentally, let us know!

    -----------
    Arg:
    X: Features/ Input vector -- torch tensor
    Y: Mapping/ Target variable vector -- torch tensor
    kernel_type: 'RBF'/'Linear'/'Periodic'/'Matern05'/'Matern15'/'Matern25' select one -- str
    noise_free: True or False (are observations noise free?)
    plot_1d: True or False (should we plot or not?)
    plot_bounds: Tuple of lower and upper bounds (xL, xU)

    ----------
    returns:

    model: a GP model object in train mode -- gpytorch
    """

    # make sure input data is shaped properly (ntrain by ninputs)
    if X.ndim != 2:
        raise ValueError(f"X must be 2D (n_train x n_inputs), got shape {X.shape}")

#     if X.shape[-1]!=2:
#         raise ValueError(f"Expected exactly 2 inputs (Polymer_thickness_mm, PCL_fraction), got {X.shape[-1]}")


    # make sure training data has the right dimension
    if Y.ndim == 1:
        Y = Y.unsqueeze(-1)
    elif Y.ndim!=2 or Y.shape[-1]!=1:
        raise ValueError(f"Y must be of shape (n,) or (n,1), got {Y.shape}")
        
    # --- determine continuous dims ---
    n_inputs = X.shape[-1]
    all_dims = list(range(n_inputs))
    cat_dims = cat_dims or []
    cont_dims = [d for d in all_dims if d not in cat_dims]
    n_cont = len(cont_dims)        
    

    # --- model transforms (internal scaling) --- # Only over continuous data
    input_transform   = Normalize(d=n_inputs, indices=cont_dims)
    outcome_transform = Standardize(m=1)
#     input_transform  = Normalize(d=2)
#     outcome_transform = Standardize(m=1)

    # --- kernel with ARD + constraints ---
    ls_constraint = Interval(ls_floor, ls_ceiling)
    ls_prior      = GammaPrior(4.0, 6.0)        # mean ≈ 0.667 (soft prior, bound dominates)
#     num_params = 2

    if kernel_type == "RBF":
        base_kernel = RBFKernel(ard_num_dims = n_cont, lengthscale_prior = ls_prior, lengthscale_constraint = ls_constraint)        
#         base_kernel = RBFKernel(ard_num_dims=num_params, lengthscale_prior=ls_prior, lengthscale_constraint=ls_constraint)
        
    elif kernel_type == "Matern25":
        base_kernel = MaternKernel(nu = 2.5, ard_num_dims = n_cont, lengthscale_prior = ls_prior, lengthscale_constraint = ls_constraint)
#         base_kernel = MaternKernel(nu=2.5, ard_num_dims=num_params, lengthscale_prior=ls_prior, lengthscale_constraint=ls_constraint)

    covar_module = ScaleKernel(base_kernel)


    # --- likelihood (noise) ---
    # likelihood (GammaPrior(shape=2.0, rate=2/mean) → mean ≈ noise_prior_mean); GammaPrior(k/r); mean ~ = k/r
    rate = 2.0 / max(noise_prior_mean, 1e-6)

    if noise_free:
        likelihood = GaussianLikelihood(
            noise_prior = GammaPrior(2.0, 40.0),   # mean ~0.05
            noise_constraint = Interval(5e-4, 0.3),
        )
    else:
        likelihood = GaussianLikelihood(
            noise_prior = GammaPrior(2.0, rate),
            noise_constraint = Interval(1e-4, 0.5),
        )

        
    # --- model ---
    
    if cat_dims:
        # cont_kernel_factory preserves our custom kernel + priors for continuous dims
        # MixedSingleTaskGP handles categorical dims with its internal CategoricalKernel
        def cont_kernel_factory(batch_shape, ard_num_dims, active_dims):
            if kernel_type == "RBF":
                base = RBFKernel(ard_num_dims = ard_num_dims, active_dims = active_dims, lengthscale_prior = ls_prior, lengthscale_constraint = ls_constraint)
            elif kernel_type == "Matern25":
                base = MaternKernel(nu = 2.5, ard_num_dims = ard_num_dims, active_dims = active_dims, lengthscale_prior = ls_prior, lengthscale_constraint = ls_constraint)
                
            return ScaleKernel(base, batch_shape=batch_shape)

        model = MixedSingleTaskGP(
            train_X = X,
            train_Y = Y,
            cat_dims = cat_dims,
            cont_kernel_factory = cont_kernel_factory,
            likelihood = likelihood,
            input_transform = input_transform,
            outcome_transform = outcome_transform
        )
    else:
        model = SingleTaskGP(
            train_X = X,
            train_Y = Y,
            covar_module = covar_module,
            likelihood = likelihood,
            input_transform = input_transform,
            outcome_transform = outcome_transform
        )        

#     model = SingleTaskGP(
#         train_X=X, # RAW physical units
#         train_Y=Y,
#         covar_module=covar_module,
#         likelihood=likelihood,
#         input_transform=input_transform,
#         outcome_transform=outcome_transform
#     )

    # --- fit ---
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(model.likelihood, model)
    from botorch.fit import fit_gpytorch_mll
    with gpytorch.settings.cholesky_jitter(1e-3):
        fit_gpytorch_mll(mll)

    model.eval()
    
#     print(model)


    # --- diagnostics ---
    if cat_dims:
        ls_add = model.covar_module.kernels[0].base_kernel.kernels[0].base_kernel.lengthscale.detach().cpu().numpy().ravel()
        ls_prod = model.covar_module.kernels[1].base_kernel.kernels[0].base_kernel.lengthscale.detach().cpu().numpy().ravel()
        print("ARD length-scales - additive term (normalized units):", ls_add)
        print("ARD length-scales - product term (normalized units):", ls_prod)
    else:
        hyp = model.covar_module.base_kernel.lengthscale.detach().cpu().numpy().ravel()
        print("ARD length-scales (normalized units):", hyp)
    print("Noise:", model.likelihood.noise.item())

    if cat_dims:
        return model, {"ARD length-scales - additive term:": ls_add, "ARD length-scales - product term:": ls_prod}, model.likelihood.noise.item()
    else:
        return model, {"ARD length-scales:": hyp}, model.likelihood.noise.item()


def run_GP_and_Plot(df, feature_cols, target_col, bounds_dict, condn_col, condn_thresh):
    
    # Detecting categorical columns and encoding them (if exists)
    cat_cols  = [c for c in feature_cols if df[c].dtype == object]
    cont_cols = [c for c in feature_cols if df[c].dtype != object]
    
    encoding_maps  = {}   # str -> int
    decoding_maps  = {}   # int -> str
    df_encoded = df.copy()
    for col in cat_cols:
        categories = sorted(df[col].unique())
        encoding_maps[col]  = {cat: i for i, cat in enumerate(categories)}
        decoding_maps[col]  = {i: cat for i, cat in enumerate(categories)}
        df_encoded[col] = df[col].map(encoding_maps[col])    
    
    # 1. Fetching input and output
    X_raw = torch.tensor(df_encoded[feature_cols].to_numpy(), dtype = torch.double)
    Y_raw = torch.tensor(df[target_col].to_numpy(), dtype = torch.double).unsqueeze(-1)   # (n,1)

    condn_val = float(df[condn_col].iloc[0])
    if condn_val == condn_thresh:
        ls_floor = 0.24      # 0.24–0.25 for “camera‑ready”
        noise_mean = 0.08
#     elif condn_val == 7.5:
#         ls_floor = 0.24      # bump to 0.24 if slight ripple remains
#         noise_mean = 0.10    # 0.10–0.12 only if ripple persists
    else:
        ls_floor = 0.20
        noise_mean = 0.08

    # 2. RAW bounds (physical space)

    cat_dims  = [feature_cols.index(c) for c in cat_cols]
    cont_dims = [feature_cols.index(c) for c in cont_cols]

    full_x_min = []
    full_x_max = []
    
    for col in feature_cols:
        if col in cat_cols:
            full_x_min.append(0.0)
            full_x_max.append(float(len(encoding_maps[col]) - 1))
        else:
            full_x_min.append(bounds_dict[col][0])
            full_x_max.append(bounds_dict[col][1])
            
    bounds_phys = torch.tensor([full_x_min, full_x_max], dtype=torch.double)    
    
#     x_min = torch.tensor([bounds_dict[c][0] for c in feature_cols], dtype = torch.double)
#     x_max = torch.tensor([bounds_dict[c][1] for c in feature_cols], dtype = torch.double)
#     bounds_phys = torch.stack([x_min, x_max])  # shape (2, 2)
    

    # 3. Train GP on raw X, Y
    model, ls_info, noise = get_trained_GP(
        X_raw, Y_raw, kernel_type = "Matern25", cat_dims = cat_dims, noise_free = False, ls_floor = ls_floor, ls_ceiling = 3.0, noise_prior_mean = noise_mean)

    # make sure model is double
    model = model.to(torch.double)

    
    # 4. Constrained qEI in RAW space

    # Build a sampler (handle older/newer BoTorch APIs)
    try:
        sampler = SobolQMCNormalSampler(sample_shape = torch.Size([1024]))
    except TypeError:
        sampler = SobolQMCNormalSampler(num_samples = 1024)

    # Objective: maximize t_cross_days in RAW units
    def obj(samples, X = None):
        return samples.squeeze(-1)  # maximize t_cross_days

    # Constraint: feasible iff t_cross_days >= 0 (in original units)
    def constraint(samples, X = None):
        return (-samples).squeeze(-1)  # feasible if t_cross_days >= 0; need negative sign because this defaults to <0 for positive sign

    constrained_obj = ConstrainedMCObjective(
        objective = obj,
        constraints = [constraint],
        infeasible_cost = torch.tensor(-1e-3, dtype = torch.double),  # small negative to keep gradient signal
    )

    feas_mask = (Y_raw.squeeze(-1) >= 0)
    best_f = Y_raw[feas_mask].max().item() if feas_mask.any() else 0.0

    ei = qExpectedImprovement(
        model = model,
        best_f = best_f,
        sampler = sampler,
        objective = constrained_obj,
    )

    # 5. Optimize EI in RAW to get next point --- mixed (continuous + categorical)

    if cat_dims:
        # Build all possible category combinations as fixed_features_list
        cat_values = [list(decoding_maps[c].keys()) for c in cat_cols]
        cat_combinations = list(itertools.product(*cat_values))
        fixed_features_list = [
            {cat_dims[i]: combo[i] for i in range(len(cat_cols))}
            for combo in cat_combinations
        ]
        candidate, acq_value = optimize_acqf_mixed(
            acq_function = ei,
            bounds = bounds_phys,
            fixed_features_list = fixed_features_list,
            q = 1,
            num_restarts = 20,
            raw_samples = 512
        )
    else:
        candidate, acq_value = optimize_acqf(
            acq_function = ei,
            bounds = bounds_phys,
            q = 1, # number of next input values to determine
            num_restarts = 20,
            raw_samples = 512
        )   
    
#     candidate, acq_value = optimize_acqf(
#         acq_function = ei,
#         bounds = bounds_phys,
#         q = 1, # number of next input values to determine
#         num_restarts = 20,
#         raw_samples = 512
#     )

    x_next = candidate[0]  # RAW

    # Converting candidate to dict
#     next_point = {feature_cols[i]: float(x_next[i].item()) for i in range(len(feature_cols))}
    
    next_point = {}
    cont_idx = 0
    for i, col in enumerate(feature_cols):
        if col in cat_cols:
            cat_int = int(round(x_next[i].item()))
            next_point[col] = decoding_maps[col][cat_int]
        else:
            next_point[col] = float(x_next[cont_idx].item())
            cont_idx += 1    
    
    print("\nSuggested next input (Unscaled):")
    for k, v in next_point.items():
        print(f"  {k} = {v:.6g}" if isinstance(v, float) else f"  {k} = {v}")
#         print(f"- {k} = {v:.6g}")
        
    print(f"\nEI at this point = {float(acq_value.item()):.6g}")

    with torch.no_grad():
        post = model.posterior(x_next.unsqueeze(0))  # shape (1,3)
        pred_mean = post.mean.item()   # model predicts time in days to cross threshold
        std = post.variance.sqrt().item()
        
    print(f"\nGP-predicted {target_col}: {pred_mean:.6g}")
    print(f"GP uncertainty (std): {std:.6g}")

    # PLOTTING
    
    vmax = round(df[df["iteration"] == 0][target_col].max() * 1.2) ### For height in plot_gp_plus_ei_2d3d
    
    figs = []
    std_ranges= []
    for i, feature in enumerate(feature_cols):
        if feature in cat_cols:
            continue  # skip categorical columns           
        fixed_values = {
            feature_cols[j]: next_point[feature_cols[j]]
            for j in range(len(feature_cols))
        }
        
        fig, std_min, std_max = plot_gp_plus_ei_slice(feature, fixed_values, model, ei, bounds_dict, df, feature_cols, target_col, condn_col)
        figs.append(fig)
        std_ranges.append((feature, std_min, std_max))
        
    fig1, fig2, fig3 = plot_gp_plus_ei_2d3d(model, ei, bounds_dict, df, target_col, condn_col, vmax, n_grid = 100)
    figs.extend([fig1, fig2, fig3])
    
    save_iteration_report(next_point, float(acq_value.item()), pred_mean, std, target_col, df, condn_col, figs, ls_info, noise, std_ranges)

    return next_point, float(acq_value.item()), pred_mean, std, figs


def plot_gp_plus_ei_slice(dim_name, fixed_vals, model, ei_acq, bounds_dict, df, feature_cols, target_col, condn_col):
    """
    dim_name : features to slice
    fixed_vals : dict of {feature:value} specifying slice location in *physical* space
    model : GP trained on *scaled* inputs and target column
    ei_acq : ExpectedImprovement object
    bounds_dict : {"feature1": (min,max), "feature2": (min,max)}
    df : dataframe with feature columns and target column
    """

    # -----------------------
    # 1) Build physical grid for the slice
    # -----------------------
    lower, upper = bounds_dict[dim_name]

    grid_phys = torch.linspace(lower, upper, 200, dtype = torch.double)

    slice_vals = {}
    for feat in feature_cols:
        if feat == dim_name:
            slice_vals[feat] = grid_phys
        else:
            slice_vals[feat] = torch.full_like(grid_phys, fixed_vals[feat], dtype = torch.double)
            
    Xtest_phys = torch.stack([slice_vals[f] for f in feature_cols], dim = 1)  # RAW
     
    # axis labeling
    x_label = dim_name
    other_feat = [f for f in feature_cols if f!= dim_name][0]
    title_value = f"{other_feat} = {fixed_vals[other_feat]:.2g}"


    # -----------------------
    # 2) GP posterior (on target_col)
    # -----------------------
    with torch.no_grad():
        post = model.posterior(Xtest_phys) # RAW → model normalizes internally
        mean = post.mean.squeeze(-1).cpu().numpy()
        std = post.variance.sqrt().squeeze(-1).cpu().numpy()
        print("\nSTD range on this slice:", std.min(), std.max())

    # upper and lower confidence bounds
    lcb = mean - 2.0 * std
    lcb = np.maximum(lcb, 0)
    ucb = mean + 2.0 * std

    # -----------------------
    # 3) EI values on this slice
    # -----------------------
    with torch.no_grad():
        ei_vals = ei_acq(Xtest_phys.unsqueeze(1)).squeeze(-1).cpu().numpy()

    # -----------------------
    # 3) Observed points projected to this axis
    # -----------------------

    x_obs = df[dim_name].to_numpy()
    y_obs = df[target_col].to_numpy()
    x_grid_np = grid_phys.cpu().numpy()

    # -----------------------
    # 6) Plot
    # -----------------------

    fig, ax1 = plt.subplots(figsize=(9, 5))

    # GP CI band
    ax1.fill_between(x_grid_np, lcb, ucb, color = "lightblue", alpha = 0.3, label = "95% CI")

    # GP mean
    ax1.plot(x_grid_np, mean, color = "blue", linewidth = 2, label = "GP mean")

    # Observed data
    observed_color, highlight_points = get_highlight_points(df, target_col)

    # ----------------------------------------------------
    # Build legend handles (same logic as main figure)
    # ----------------------------------------------------
    legend_handles = []

    # Observed points (LHS)
    h_obs = ax1.scatter([], [], c = observed_color, marker = "x", s = 40, label = "LHS")
    legend_handles.append(h_obs)

    # Highlight points (LHS best, GP design1, GP design2, ...)
    iteration_max = df["iteration"].max()
    for i, (df_high, color) in enumerate(highlight_points):
        if iteration_max == 0:
            label = "LHS best"
        else:
            label = "LHS best" if i == 0 else f"GP design{i}"
        h = ax1.scatter([], [], c = color, marker = "x", s = 90, linewidths = 2.5, label = label)
        
        legend_handles.append(h)

    # Observed points (mirrored color)
    ax1.scatter(
        x_obs, y_obs,
        c = observed_color, s = 40, marker = "x")

    # Highlight points (mirrored colors)
    for df_high, color in highlight_points:
        ax1.scatter(
            df_high[dim_name], df_high[target_col],
            c = color, s = 90, marker = "x", linewidths = 2.5)

    ax1.set_xlabel(x_label)
    ax1.set_ylabel(target_col)
    ax1.grid(True)

    # EI on twin axis
    ax2 = ax1.twinx()
    ax2.plot(x_grid_np, ei_vals, color = "orange", linestyle = "--", linewidth = 2, label = "EI")
    ax2.set_ylabel("EI", color = "orange")

    # ----------------------------------------------------
    # Shared bottom legend (2 rows, full width)
    # ----------------------------------------------------
    n_handles = len(legend_handles)
    ncol = int(np.ceil(n_handles / 3))   # split into 3 rows

    fig.legend(
        handles = legend_handles,
        loc = "lower center",
        bbox_to_anchor = (0.5, -0.18),
        ncol = ncol,
        frameon = True
    )

    plt.subplots_adjust(bottom = 0.28)   # reserve space for 2-row legend

    plt.title(f"GP Posterior + EI Slice at {title_value}")
    fig.canvas.draw()
    plt.tight_layout()
    ax2.set_ylim(0, ei_vals.max()*1.05)
    # align zeroes between ax1 and ax2
    ei_cap = float(np.nanmax(ei_vals)) * 1.05 if np.isfinite(np.nanmax(ei_vals)) else None
    align_yzeros_pixel_exact(ax1, ax2, y_left = 0.0, y_right = 0.0,
                            ensure_zero_in_left = True,
                            post_pad_right_max = ei_cap)

    save_plot(fig, f"GP_EI_slice_{dim_name}", df, condn_col)
    plt.show()
    
    return fig, std.min(), std.max()

    
    
def plot_gp_plus_ei_2d3d(model, ei_acq, bounds_dict, df, target_col, condn_col, vmax, n_grid = 100):
    """
    model: GP trained on scaled inputs, predicting target_col
    ei_acq: ExpectedImprovement object
    bounds_dict: {"feature1": (min,max), "feature2": (min,max)}
    df: dataframe with feature columns and target column
    x_min, x_max: tensors used for scaling
    """

    feature_cols = list(bounds_dict.keys())
    x_feat = feature_cols[1]
    y_feat = feature_cols[0]
    
    # -----------------------
    # 1) Build 2D physical grid
    # -----------------------
    th_min, th_max = bounds_dict[y_feat]
    pcl_min, pcl_max = bounds_dict[x_feat]

    th_grid = torch.linspace(th_min, th_max, n_grid, dtype = torch.double)
    pcl_grid = torch.linspace(pcl_min, pcl_max, n_grid, dtype = torch.double)

    TH, PCL = torch.meshgrid(th_grid, pcl_grid, indexing = "ij")  # (n,n)

    X_phys = torch.stack([TH.reshape(-1), PCL.reshape(-1)], dim = 1)  # RAW

    # -----------------------
    # 2) GP posterior & EI valuespl
    # -----------------------
    with torch.no_grad():
        post = model.posterior(X_phys)
        mean = post.mean.view(-1).cpu().numpy()
        std = post.variance.sqrt().view(-1).cpu().numpy()
        ei_vals = ei_acq(X_phys.unsqueeze(1)).view(-1).cpu().numpy()

    # -----------------------
    # 3) Reshape for plotting
    # -----------------------
    mean_2d = mean.reshape(n_grid, n_grid)
    std_2d = std.reshape(n_grid, n_grid)
    ei_2d = ei_vals.reshape(n_grid, n_grid)

    # -----------------------
    # 4) 2D mean and EI
    # -----------------------
    fig1, ax = plt.subplots(1, 2, figsize=(14, 6))
    # After computing mean_2d for each dataset, gather global limits:
    #vmin = min(mean_2d_5.min(),  mean_2d_7pt5.min())
    #vmax = max(mean_2d_5.max(),  mean_2d_7pt5.max())
    vmin = 0
#     vmax = 450
    im0 = ax[0].imshow(
        mean_2d,
        origin = "lower",
        extent = [pcl_min, pcl_max, th_min, th_max],
        aspect = "auto",
        cmap = "viridis",
        vmin = vmin, vmax = vmax
    )
    ax[0].set_title(f"GP Posterior Mean ({target_col})")
    ax[0].set_xlabel(x_feat)
    ax[0].set_ylabel(y_feat)

    fig1.colorbar(im0, ax = ax[0])

    # ----------------------------------------------------
    # Plot all observed points (red)
    # ----------------------------------------------------
    ax[0].scatter(
        df[x_feat], df[y_feat],
        c = "red", s = 40, marker = "x", label = "LHS")

    cs = ax[0].contour(np.linspace(pcl_min, pcl_max, n_grid),
                   np.linspace(th_min, th_max, n_grid),
                   mean_2d, levels = 10, colors = "k",
                   linewidths = 0.6, alpha = 0.5)
    
    ax[0].clabel(cs, inline = True, fmt = "%.0f", fontsize = 9)

    # ----------------------------------------------------
    # NEW LOGIC: highlight based on target column
    # ----------------------------------------------------
    iteration_max = df["iteration"].max()

    # A palette with enough distinct colors
    color_list = ["cyan", "magenta", "lime", "orange", "blue", "brown", "silver", "gold", "green", "chocolate", "orchid"]

    # We will collect legend handles manually
    legend_handles = []

    if iteration_max == 0:
        # Only one iteration → highlight the single highest target_col
        idx = df[target_col].idxmax()
        df_high = df.loc[[idx]]

        h = ax[0].scatter(
            df_high[x_feat], df_high[y_feat],
            c = "cyan", s = 90, marker = "x", linewidths = 2.5, label = "LHS best")
        
        legend_handles.append(h)

    else:
        # Multiple iterations → highlight the max target_col within each iteration
        for ds_id in range(0, iteration_max + 1):
            df_sub = df[df["iteration"] == ds_id]
            if len(df_sub) == 0:
                continue

            idx = df_sub[target_col].idxmax()
            df_high = df.loc[[idx]]

            color = color_list[ds_id % len(color_list)]

            # Labeling rule:
            if ds_id == 0:
                label = "LHS best"
            else:
                label = f"GP design{ds_id}"

            h = ax[0].scatter(
                df_high[x_feat], df_high[y_feat],
                c = color, s = 90, marker = "x", linewidths = 2.5, label = label)
            
            legend_handles.append(h)


    # -----------------------
    # 7) Plot EI
    # -----------------------
    im1 = ax[1].imshow(ei_2d, origin = "lower",
                       extent = [pcl_min, pcl_max, th_min, th_max],
                       aspect = "auto", cmap = "inferno")
    ax[1].set_title("Expected Improvement")
    ax[1].set_xlabel(x_feat)
    ax[1].set_ylabel(y_feat)
    fig1.colorbar(im1, ax = ax[1])

    # Overlay observed points
    ax[1].scatter(df[x_feat], df[y_feat],
                  c = "red", s = 40, marker = "x")

    # ----------------------------------------------------
    # Mirror highlight colors on EI panel
    # ----------------------------------------------------
    if iteration_max == 0:
        # Single iteration → highlight the single highest target
        idx = df[target_col].idxmax()
        df_high = df.loc[[idx]]

        ax[1].scatter(
            df_high[x_feat], df_high[y_feat],
            c = "cyan", # same highlight color as left panel
            s = 90, marker = "x", linewidths=2.5)

    else:
        # Multiple iterations → highlight max target within each iteration
        for ds_id in range(0, iteration_max + 1):
            df_sub = df[df["iteration"] == ds_id]
            if len(df_sub) == 0:
                continue

            idx = df_sub[target_col].idxmax()
            df_high = df.loc[[idx]]

            color = color_list[ds_id % len(color_list)]

            ax[1].scatter(
                df_high[x_feat], df_high[y_feat],
                c = color, # same per‑iteration color as left panel
                s = 90, marker = "x", linewidths = 2.5)


    # ----------------------------------------------------
    # Shared legend for both subplots (bottom, 2 rows)
    # ----------------------------------------------------
    all_handles = [ax[0].collections[0]] + legend_handles
    n_handles = len(all_handles)
    ncol = int(np.ceil(n_handles / 3))   # split into 3 rows

    plt.tight_layout()

    fig1.legend(
        handles = all_handles,
        loc = "lower center",
        bbox_to_anchor = (0.5, -0.10),   # slightly lower to fit 2 rows
        ncol = ncol,
        frameon = True
    )

    plt.subplots_adjust(bottom = 0.22)   # reserve more space for 2-row legend
    save_plot(fig1, "GP_2D_mean_EI", df, condn_col)

    """
    TH, PCL: meshgrid tensors (n,n)
    mean_2d, std_2d: numpy arrays (n,n)
    """

    # 95% confidence interval surfaces
    lcb_2d = mean_2d - 2 * std_2d
    lcb_2d = np.maximum(lcb_2d, 0)
    ucb_2d = mean_2d + 2 * std_2d



    # Convert to numpy for plotting
    TH_np = TH.cpu().numpy()
    PCL_np = PCL.cpu().numpy()

    fig2 = plt.figure(figsize = (12, 8))
    ax = fig2.add_subplot(111, projection = "3d")
    #ax.set_proj_type('ortho')
    # explicit limits & aspect (prevents autoscale surprises)
    # After computing mean_2d for each dataset, gather global limits:
    #zmin = min(mean_2d_5.min(),  mean_2d_7pt5.min())
    #zmax = max(mean_2d_5.max(),  mean_2d_7pt5.max())
    zmin = 0
    zmax = vmax
    ax.set_xlim(pcl_min, pcl_max)
    ax.set_ylim(th_min, th_max)
    ax.set_zlim(zmin, zmax)


    # Mean surface
    surf = ax.plot_surface(PCL_np, TH_np, mean_2d,
                    cmap = "viridis", alpha = 0.9, linewidth = 0, antialiased = True)

    # Lower CI surface
    #ax.plot_surface(PCL_np, TH_np, lcb_2d,
    #                color="blue", alpha=0.25, linewidth=0)

    # Upper CI surface
    #ax.plot_surface(PCL_np, TH_np, ucb_2d,
    #                color="blue", alpha=0.25, linewidth=0)

    # Put colorbar in its own inset axes so it doesn't squeeze the 3D plot
    cax = fig2.add_axes([0.825, 0.20, 0.02, 0.60])  # [left, bottom, width, height] in figure coords
    surf.set_clim(vmin = zmin, vmax = zmax)
    cb = fig2.colorbar(surf, cax = cax)
    cb.ax.tick_params(labelsize = 10)


    # Optionally set a nice view angle
    #ax.view_init(elev=25, azim=-60)


    ax.set_xlabel(x_feat,labelpad = 10, fontsize = 12)
    ax.set_ylabel(y_feat,labelpad = 10, fontsize = 12)
    ax.set_zlabel(target_col, labelpad = 10, fontsize = 12)
    #ax.set_title("3D GP Posterior with 95% Confidence Surfaces")
    #ax.set_title("3D GP Posterior")



    plt.show()
    save_plot(fig2, "GP_3D_posterior", df, condn_col)



    """
    TH, PCL: meshgrid tensors (n,n)
    ei_2d: numpy array (n,n) of EI values
    """

    TH_np = TH.cpu().numpy()
    PCL_np = PCL.cpu().numpy()

    fig3 = plt.figure(figsize = (12, 8))
    ax = fig3.add_subplot(111, projection = "3d")
    #ax.set_proj_type('ortho')

    ax.set_xlim(pcl_min, pcl_max)
    ax.set_ylim(th_min, th_max)

    ax.set_zlim(0, float(np.nanmax(ei_2d)) * 1.05)

    #ax.set_box_aspect((1,1,0.7))  # keep same 0.6 as above



    # EI surface
    surf = ax.plot_surface(PCL_np, TH_np, ei_2d,
                          cmap = "inferno", linewidth = 0,
                          antialiased = True, alpha = 0.9)

    # Axis labels with padding and font size
    ax.set_xlabel(x_feat, labelpad = 10, fontsize = 12)
    ax.set_ylabel(y_feat, labelpad = 10, fontsize = 12)
    ax.set_zlabel("Expected Improvement", labelpad = 10, fontsize = 12)  # Increased labelpad

    # Title
    #ax.set_title("3D Expected Improvement Surface", pad=20)


    cax = fig3.add_axes([0.825, 0.20, 0.02, 0.60])
    cb = fig3.colorbar(surf, cax = cax)
    cb.ax.tick_params(labelsize = 10)

    #ax.view_init(elev=25, azim=-60)


    plt.show()
    save_plot(fig3, "GP_3D_EI_surface", df, condn_col)
    
    return fig1, fig2, fig3

