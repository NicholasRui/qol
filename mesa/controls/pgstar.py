# Methods handling controls about pgstar

def show_pgstar(self, abs_path=None):
    self.add_control(namelist='star_job', category='enable pgstar',
            control='pgstar_flag', value=True)

def pgstar_Grid1_enable(self, num_cols, num_rows, reset=False,
        Grid1_win_width=14, Grid1_win_aspect_ratio=0.5,
        write_path=None,
        Grid1_file_prefix='Grid1_', Grid1_file_interval=5, Grid1_file_width=-1, Grid1_file_aspect_ratio=-1,
        ):
    """
    Enable use of Grid1 pgstar plot

    if writedir is not None, save plots to that path

    Make animation with:
    ffmpeg -r 10 -pattern_type glob -i "Grid1_*.png" -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -profile:v high -crf 25 -pix_fmt yuv420p Grid1_anim.mp4
    """
    self.add_control(namelist='pgstar', category='Grid1',
            control='Grid1_win_flag', value=True)
    self.add_control(namelist='pgstar', category='Grid1',
            control='Grid1_win_width', value=Grid1_win_width)
    self.add_control(namelist='pgstar', category='Grid1',
            control='Grid1_win_aspect_ratio', value=Grid1_win_aspect_ratio)
    self.add_control(namelist='pgstar', category='Grid1',
            control='Grid1_num_cols', value=num_cols)
    self.add_control(namelist='pgstar', category='Grid1',
            control='Grid1_num_rows', value=num_rows)
    self.add_control(namelist='pgstar', category='Grid1', comment='<= 10',
            control='Grid1_num_plots', value=num_cols*num_rows)

    if write_path is not None:
        if write_path[-1] != '/':
           write_path += '/'
        
        self.add_control(namelist='pgstar', category='Grid1: save',
                control='Grid1_file_flag', value=True)
        self.add_control(namelist='pgstar', category='Grid1: save',
                control='Grid1_file_dir', value=write_path)
        self.add_control(namelist='pgstar', category='Grid1: save',
                control='Grid1_file_prefix', value=Grid1_file_prefix)
        self.add_control(namelist='pgstar', category='Grid1: save',
                control='Grid1_file_interval', value=Grid1_file_interval) # ! output when mod(model_number,Grid1_file_interval)==0
        self.add_control(namelist='pgstar', category='Grid1: save',
                control='Grid1_file_width', value=Grid1_file_width) # (inches) negative means use same value as for window
        self.add_control(namelist='pgstar', category='Grid1: save',
                control='Grid1_file_aspect_ratio', value=Grid1_file_aspect_ratio)
    
    # Reset defaults
    if reset:
        # note: the category actually matters here: this reset needs to be sorted before any new plots,
        #       otherwise it erases everything
        self.add_control(namelist='pgstar', category='Grid1',
                control='Grid1_plot_name(:)', value='')
        self.add_control(namelist='pgstar', category='Grid1', comment='number from 1 at top',
                control='Grid1_plot_row(:)', value=1)
        self.add_control(namelist='pgstar', category='Grid1', comment='plot spans this number of rows',
                control='Grid1_plot_rowspan(:)', value=1)
        self.add_control(namelist='pgstar', category='Grid1', comment='number from 1 at left',
                control='Grid1_plot_col(:)', value=1)
        self.add_control(namelist='pgstar', category='Grid1', comment='plot spans this number of columns',
                control='Grid1_plot_colspan(:)', value=1)

def pgstar_Grid1_add_plot(self, plot_col, plot_row, plot_name,
            pad_left=0., pad_right=0., pad_top=0., pad_bot=0., txt_scale_factor=1.
            ):
    """
    Add a panel to Grid1, making sure that the position does not already exist
    """
    assert (plot_col, plot_row) not in self.Grid1_plot_coords
    self.Grid1_plot_coords.append((plot_col, plot_row))

    plot_num = len(self.Grid1_plot_coords)

    category = f'Grid1({plot_num}): {plot_name}'
    self.add_control(namelist='pgstar', category=category,
            control=f'Grid1_plot_name({plot_num})', value=plot_name)
    self.add_control(namelist='pgstar', category=category,
            control=f'Grid1_plot_row({plot_num})', value=plot_row)
    self.add_control(namelist='pgstar', category=category,
            control=f'Grid1_plot_col({plot_num})', value=plot_col)
    self.add_control(namelist='pgstar', category=category,
            control=f'Grid1_plot_pad_left({plot_num})', value=pad_left)
    self.add_control(namelist='pgstar', category=category,
            control=f'Grid1_plot_pad_right({plot_num})', value=pad_right)
    self.add_control(namelist='pgstar', category=category,
            control=f'Grid1_plot_pad_top({plot_num})', value=pad_top)
    self.add_control(namelist='pgstar', category=category,
            control=f'Grid1_plot_pad_bot({plot_num})', value=pad_bot)
    self.add_control(namelist='pgstar', category=category,
            control=f'Grid1_txt_scale_factor({plot_num})', value=txt_scale_factor)

    # Assume single panel
    self.add_control(namelist='pgstar', category=category,
            control=f'Grid1_plot_colspan({plot_num})', value=1)
    self.add_control(namelist='pgstar', category=category,
            control=f'Grid1_plot_rowspan({plot_num})', value=1)

def pgstar_Grid1_qol_default(self, write_path):
    """
    Implement default qol settings for Grid1
    """
    self.show_pgstar()

    self.pgstar_Grid1_enable(num_cols=3, num_rows=2,
        Grid1_win_width=14, Grid1_win_aspect_ratio=0.5,
        write_path=write_path,
        Grid1_file_prefix='Grid1_', Grid1_file_interval=5, Grid1_file_width=-1, Grid1_file_aspect_ratio=-1)

    self.pgstar_Grid1_add_plot(plot_col=1, plot_row=1, plot_name='HR',
            pad_left=0., pad_right=0.06, pad_top=0., pad_bot=0.06, txt_scale_factor=0.7)
    self.pgstar_Grid1_add_plot(plot_col=1, plot_row=2, plot_name='TRho_Profile',
            pad_left=0., pad_right=0.06, pad_top=0.06, pad_bot=0., txt_scale_factor=0.7)
    self.pgstar_Grid1_add_plot(plot_col=2, plot_row=1, plot_name='Abundance',
            pad_left=0.02, pad_right=0.06, pad_top=0., pad_bot=0.06, txt_scale_factor=0.7)
    self.pgstar_Grid1_add_plot(plot_col=2, plot_row=2, plot_name='Power',
            pad_left=0.02, pad_right=0.06, pad_top=0.06, pad_bot=0., txt_scale_factor=0.7)

    category = 'Grid1(5): Profile_Panels1'
    self.pgstar_Grid1_add_plot(plot_col=3, plot_row=1, plot_name='Profile_Panels1',
            pad_left=0.06, pad_right=0.02, pad_top=0., pad_bot=0.06, txt_scale_factor=0.7)
    self.add_control(namelist='pgstar', category=category,
            control=f'Profile_Panels1_num_panels', value=1)
    self.add_control(namelist='pgstar', category=category,
            control=f'Profile_Panels1_yaxis_name(1)', value='velocity')
    self.add_control(namelist='pgstar', category=category,
            control=f'Profile_Panels1_other_yaxis_name(1)', value='entropy')

    category = f'Grid1(6): Text_Summary1'
    self.pgstar_Grid1_add_plot(plot_col=3, plot_row=2, plot_name='Text_Summary1',
            pad_left=0.06, pad_right=0., pad_top=0.06, pad_bot=0., txt_scale_factor=0.3)
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_num_rows', value=12)
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_num_cols', value=1)
    self.add_control(namelist='pgstar', category=category, comment='reset to default',
            control='Text_Summary1_name(:,:)', value='')
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(1,1)', value='star_mass')
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(2,1)', value='he_core_mass')
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(3,1)', value='co_core_mass')
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(4,1)', value='log_L')
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(5,1)', value='log_LH')
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(6,1)', value='log_LHe')
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(7,1)', value='log_LZ')    
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(8,1)', value='log_Teff')
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(9,1)', value='log_R')
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(10,1)', value='log_star_age')
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(11,1)', value='log_dt')
    self.add_control(namelist='pgstar', category=category,
            control='Text_Summary1_name(12,1)', value='num_zones')

